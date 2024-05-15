classdef stewart < nosnoc.dcs.base
    properties
        model

        theta
        theta_sys
        lambda
        lambda_sys
        mu
        mu_sys

        z_all

        f_x
        g_Stewart

        sys_idx

        dims

        % TODO some of these should be in base perhaps
        f_x_fun
        f_q_fun
        g_z_fun
        g_alg_fun
        g_switching_fun
        g_Stewart_fun
        lambda00_fun
        g_path_fun
        g_comp_path_fun
        g_terminal_fun
        f_q_T_fun
        f_lsq_x_fun
        f_lsq_u_fun
        f_lsq_T_fun
    end

    methods
        function obj = stewart(model)
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

            % symbolic variables z = [theta;lambda;mu_Stewart];
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
                g_switching = [g_switching; model.g_ind{ii} - obj.lambda_sys{ii}+obj.mu_sys{ii}*ones(dims.n_f_sys(ii),1)];
                g_convex = [g_convex;ones(dims.n_f_sys(ii),1)'*obj.theta_sys{ii} - 1];
                lambda00_expr = [lambda00_expr; model.g_ind{ii} - min(model.g_ind{ii})];
            end
            g_alg = [g_switching;g_convex];

            obj.f_x_fun = Function('f_x', {model.x, model.z, obj.lambda, obj.theta, obj.mu, model.u, model.v_global, model.p}, {obj.f_x, model.f_q});
            obj.f_q_fun = Function('f_q', {model.x, model.z, obj.lambda, obj.theta, obj.mu, model.u, model.v_global, model.p}, {model.f_q});
            obj.g_z_fun = Function('g_z', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_z});
            obj.g_alg_fun = Function('g_alg', {model.x, model.z, obj.lambda, obj.theta, obj.mu, model.u, model.v_global, model.p}, {g_alg});
            obj.g_switching_fun = Function('g_switching', {model.x, model.z, obj.lambda, obj.mu, model.v_global, model.p}, {g_switching});
            obj.g_Stewart_fun = Function('g_Stewart', {model.x, model.z, model.v_global, model.p}, {model.g_ind{:}});
            obj.lambda00_fun = Function('lambda00', {model.x, model.z, model.v_global, model.p_global}, {lambda00_expr});
            obj.g_path_fun = Function('g_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_path}); % TODO(@anton) do dependence checking for spliting the path constriants
            obj.g_comp_path_fun  = Function('g_comp_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_comp_path});
            obj.g_terminal_fun  = Function('g_terminal', {model.x, model.z, model.v_global, model.p_global}, {model.g_terminal});
            obj.f_q_T_fun = Function('f_q_T', {model.x, model.z, model.v_global, model.p}, {model.f_q_T});
            obj.f_lsq_x_fun = Function('f_lsq_x_fun',{model.x,model.x_ref,model.p},{model.f_lsq_x});
            obj.f_lsq_u_fun = Function('f_lsq_u_fun',{model.u,model.u_ref,model.p},{model.f_lsq_u});
            obj.f_lsq_T_fun = Function('f_lsq_T_fun',{model.x,model.x_ref_end,model.p_global},{model.f_lsq_T});
        end
    end
end
