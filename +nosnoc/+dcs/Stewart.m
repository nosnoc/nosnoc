classdef Stewart < nosnoc.dcs.Base
    properties

        theta % casadi.SX|casadi.MX: Convex multipliers for the Stewart LP.
        theta_sys % cell(casadi.SX|casadi.MX): Cell array with each cell containing theta for linearly independent subsystem.
        lambda % casadi.SX|casadi.MX: Inequality lagrange multiplier in Stewart LP.
        lambda_sys  % cell(casadi.SX|casadi.MX): Cell array with each cell containing lambda for linearly independent subsystem.
        mu % casadi.SX|casadi.MX: Equality lagrange multiplier in Stewart LP.
        mu_sys % cell(casadi.SX|casadi.MX): Cell array with each cell containing mu for linearly independent subsystem.

        z_all % casadi.SX|casadi.MX: Vector of all algorithmic algebraic variables.

        f_x  % casadi.SX|casadi.MX: r.h.s. of the ODE, f_x = sum_i F_i*theta_i , i is the index of the subystems.
        g_indicator % casadi.SX|casadi.MX: Expression for Stewart indicator functions g(x).
        g_alg % casadi.SX|casadi.MX: Expression for algebraic part of DCS.

        dims % struct: Contains relevant dimensions.
        
        
        g_lp_stationarity_fun % casadi.Function: Stewart LP stationarity function.
        g_indicator_fun % casadi.Function: Function for Stewart indicator functions.
        lambda00_fun % casadi.Function: Function mapping initial lambda from initial x.
    end

    methods
        function obj = Stewart(model)
            obj.model = model;
            obj.dims = model.dims;
        end
        
        function generate_variables(obj, opts)
            import casadi.*
            dims = obj.dims;
            % dimensions
            dims.n_theta = sum(obj.dims.n_f_sys); % number of modes
            dims.n_lambda = dims.n_theta;
            dims.n_mu = dims.n_sys;
            idx = 1;
            for ii = 1:dims.n_sys
                sys_idx{ii} = idx:(idx + obj.dims.n_f_sys(ii)-1);
                idx = idx + obj.dims.n_f_sys(ii);
                ii_str = num2str(ii);
                % define theta (Filippov multiplers)
                obj.theta_sys{ii} = define_casadi_symbolic(opts.casadi_symbolic_mode,['theta_' ii_str],obj.dims.n_f_sys(ii));
                obj.theta = [obj.theta;obj.theta_sys{ii}];
                % define mu_i (Lagrange multipler of e'theta =1;)
                obj.mu_sys{ii} = define_casadi_symbolic(opts.casadi_symbolic_mode,['mu_' ii_str],1);
                obj.mu = [obj.mu;obj.mu_sys{ii}];
                % define lambda_i (Lagrange multipler of theta >= 0;)
                obj.lambda_sys{ii} = define_casadi_symbolic(opts.casadi_symbolic_mode,['lambda_' ii_str],obj.dims.n_f_sys(ii));
                obj.lambda = [obj.lambda;obj.lambda_sys{ii}];
            end

            % symbolic variables z_all = [theta;lambda_stewart;mu_stewart;z_user_algebarics];
            obj.z_all = [obj.theta;obj.lambda;obj.mu;obj.model.z];
            obj.dims = dims;
        end

        function generate_equations(obj, opts)
            import casadi.*
            model = obj.model;
            dims = obj.dims;

            
            obj.f_x = zeros(dims.n_x,1);
            for ii = 1:dims.n_sys
                obj.f_x = obj.f_x + model.F{ii}*obj.theta_sys{ii};
            end
            obj.g_indicator = model.g_indicator;

            g_lp_stationarity = []; % collects lp_stationarity function algebraic equations, 0 = g_i(x) - \lambda_i - e \mu_i
            g_convex = []; % equation for the convex multiplers 1 = e' \theta
            
            lambda00_expr =[];
            for ii = 1:dims.n_sys
                % algebraic equations and complementarity condtions of the DCS
                % (Note that the cross complementarities are later defined when the discrete
                % time variables for every IRK stage in the create_nlp_nosnoc function are defined.)
                % g_indicator_i - lambda_i + mu_i e_i = 0; for all i = 1,..., n_sys
                % lambda_i'*theta_i = 0; for all i = 1,..., n_sys
                % lambda_i >= 0;    for all i = 1,..., n_sys
                % theta_i >= 0;     for all i = 1,..., n_sys
                % Gradient of Lagrange Function of indicator LP
                g_lp_stationarity = [g_lp_stationarity; model.g_indicator{ii} - obj.lambda_sys{ii}+obj.mu_sys{ii}*ones(dims.n_f_sys(ii),1)];
                g_convex = [g_convex;ones(dims.n_f_sys(ii),1)'*obj.theta_sys{ii} - 1];
                lambda00_expr = [lambda00_expr; model.g_indicator{ii} - min(model.g_ind{ii})];
            end
            g_alg = [g_lp_stationarity;g_convex];
            obj.g_alg = g_alg;

            obj.f_x_fun = Function('f_x', {model.x, model.z, obj.lambda, obj.theta, obj.mu, model.u, model.v_global, model.p}, {obj.f_x, model.f_q});
            obj.f_q_fun = Function('f_q', {model.x, model.z, obj.lambda, obj.theta, obj.mu, model.u, model.v_global, model.p}, {model.f_q});
            obj.g_z_fun = Function('g_z', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_z});
            obj.g_alg_fun = Function('g_alg', {model.x, model.z, obj.lambda, obj.theta, obj.mu, model.u, model.v_global, model.p}, {g_alg});
            obj.g_lp_stationarity_fun = Function('g_lp_stationarity', {model.x, model.z, obj.lambda, obj.mu, model.v_global, model.p}, {g_lp_stationarity});
            obj.g_indicator_fun = Function('g_indicator', {model.x, model.z, model.v_global, model.p}, {model.g_indicator{:}});
            obj.lambda00_fun = Function('lambda00', {model.x, model.z, model.v_global, model.p_global}, {lambda00_expr});
            obj.g_path_fun = Function('g_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_path}); % TODO(@anton) do dependence checking for spliting the path constriants
            obj.G_path_fun  = Function('G_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.G_path});
            obj.H_path_fun  = Function('H_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.H_path});
            obj.g_terminal_fun  = Function('g_terminal', {model.x, model.z, model.v_global, model.p_global}, {model.g_terminal});
            obj.f_q_T_fun = Function('f_q_T', {model.x, model.z, model.v_global, model.p}, {model.f_q_T});
            obj.f_lsq_x_fun = Function('f_lsq_x_fun',{model.x,model.x_ref,model.p},{model.f_lsq_x});
            obj.f_lsq_u_fun = Function('f_lsq_u_fun',{model.u,model.u_ref,model.p},{model.f_lsq_u});
            obj.f_lsq_T_fun = Function('f_lsq_T_fun',{model.x,model.x_ref_end,model.p_global},{model.f_lsq_T});
        end
    end

    methods(Access=protected)
        function propgrp = getPropertyGroups(obj)
            propgrp = getPropertyGroups@nosnoc.dcs.Base(obj);
            group_title = 'Variables';
            var_list = struct;
            var_list.x = obj.model.x;
            if ~any(size(obj.model.u) == 0)
                var_list.u = obj.model.u;
            end
            if ~any(size(obj.model.z) == 0)
                var_list.z = obj.model.z;
            end
            var_list.lambda = obj.lambda;
            var_list.theta = obj.theta;
            var_list.mu = obj.mu;
            propgrp(2) = matlab.mixin.util.PropertyGroup(var_list, group_title);
        end
    end
end
