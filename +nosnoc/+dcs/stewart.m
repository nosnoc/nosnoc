classdef stewart < nosnoc.dcs.base
    properties
        model nosnoc.model.pss

        theta
        theta_sys
        lambda
        lambda_sys
        mu
        mu_sys

        f_x
        g_Stewart

        dims
    end

    methods
        function obj = stewart(model)
            obj.model = model;
            dims = model.dims;
        end
        
        function generate_variables(obj, opts)
            import casadi.*
            dims = obj.dims;
            % dimensions
            dims.n_theta = sum(obj.dims.n_f_sys); % number of modes
            dims.n_lambda = dims.n_theta;
            for ii = 1:dims.n_sys
                ii_str = num2str(ii);
                % define theta (Filippov multiplers)
                obj.theta_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['theta_' ii_str],obj.dims.n_f_sys(ii));
                obj.theta = [obj.theta;obj.theta_sys{ii}];
                % define mu_i (Lagrange multipler of e'theta =1;)
                obj.mu_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['mu_' ii_str],1);
                obj.mu = [obj.mu;obj.mu_sys{ii}];
                % define lambda_i (Lagrange multipler of theta >= 0;)
                obj.lambda_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['lambda_' ii_str],obj.dims.n_f_sys(ii));
                obj.lambda = [obj.lambda;obj.lambda_sys{ii}];
            end

            % symbolic variables z = [theta;lambda;mu_Stewart];
            obj.z_all = [obj.theta;obj.lambda;obj.mu;obj.z];
        end

        function generate_equations(obj, opts)
            import casadi.*
            model = obj.model
            dims = obj.dims;
            
            obj.f_x = zeros(dims.n_x,1);
            for ii = 1:dims.n_sys
                obj.f_x = obj.f_x + obj.F{ii}*obj.theta_sys{ii};
                obj.g_Stewart{ii} = -obj.S{ii}*obj.c{ii};
            end

            g_switching = []; % collects switching function algebraic equations, 0 = g_i(x) - \lambda_i - e \mu_i, 0 = c(x)-lambda_p+lambda_n
            g_convex = []; % equation for the convex multiplers 1 = e' \theta
            
            lambda00_expr =[];
            for ii = 1:dims.n_sys
                % basic algebraic equations and complementarity condtions of the DCS
                % (Note that the cross complementarities are later defined when the discrete
                % time variables for every IRK stage in the create_nlp_nosnoc function are defined.)
                % g_ind_i - lambda_i + mu_i e_i = 0; for all i = 1,..., n_sys
                % lambda_i'*theta_i = 0; for all i = 1,..., n_sys
                % lambda_i >= 0;    for all i = 1,..., n_sys
                % theta_i >= 0;     for all i = 1,..., n_sys
                % Gradient of Lagrange Function of indicator LP
                g_switching = [g_switching; obj.g_Stewart{ii} - obj.lambda_sys{ii}+obj.mu_sys{ii}*ones(dims.n_f_sys(ii),1)];
                g_convex = [g_convex;ones(dims.n_f_sys(ii),1)'*obj.theta_sys{ii} - 1];
                lambda00_expr = [lambda00_expr; obj.g_Stewart{ii} - min(obj.g_Stewart{ii})];
            end
            g_alg = [g_switching;g_convex];

            obj.f_x_fun = Function('f_x', {obj.x, obj.z, obj.lambda, obj.theta, obj.mu, obj.u, obj.v_global, obj.p}, {obj.f_x, model.f_q});
            obj.f_q_fun = Function('f_q', {obj.x, obj.z, obj.lambda, obj.theta, obj.mu, obj.u, obj.v_global, obj.p}, {model.f_q});
            obj.g_z_fun = Function('g_z', {obj.x, obj.z, obj.u, obj.v_global, obj.p}, {model.g_z});
            obj.g_alg_fun = Function('g_alg', {obj.x, obj.z, obj.lambda, obj.theta, obj.mu, obj.u, obj.v_global, obj.p}, {g_alg});
            obj.g_Stewart_fun = Function('g_Stewart', {obj.x, obj.z, obj.v_global, obj.p}, {obj.g_Stewart{:}});
            obj.lambda00_fun = Function('lambda00', {obj.x, obj.z, obj.v_global, obj.p_global}, {lambda00_expr});
            obj.g_path_fun = Function('g_path', {obj.x, obj.z, obj.u, obj.v_global, obj.p}, {model.g_path}); % TODO(@anton) do dependence checking for spliting the path constriants
            obj.g_comp_path_fun  = Function('g_comp_path', {obj.x, obj.z, obj.u, obj.v_global, obj.p}, {obj.g_comp_path});
            obj.g_terminal_fun  = Function('g_terminal', {obj.x, obj.z, obj.v_global, obj.p_global}, {obj.g_terminal});
            obj.f_q_T_fun = Function('f_q_T', {obj.x, obj.z, obj.v_global, obj.p}, {obj.f_q_T});
            obj.f_lsq_x_fun = Function('f_lsq_x_fun',{obj.x,obj.x_ref,obj.p},{obj.f_lsq_x});
            obj.f_lsq_u_fun = Function('f_lsq_u_fun',{obj.u,obj.u_ref,obj.p},{obj.f_lsq_u});
            obj.f_lsq_T_fun = Function('f_lsq_T_fun',{obj.x,obj.x_ref_end,obj.p_global},{obj.f_lsq_T});
        end
    end
end
