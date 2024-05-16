classdef Heaviside < nosnoc.dcs.Base
    properties
        model

        alpha
        alpha_sys
        lambda_n
        lambda_n_sys
        lambda_p
        lambda_p_sys
        % These are relevant only for lifting. For now wait to implement until re-implementing time freezing
        % beta 
        % gamma
        % theta
        % theta_sys
        theta_expr_sys

        z_all

        f_x
        g_Stewart

        sys_idx

        dims

        g_switching_fun
        lambda00_fun
    end

    methods
        function obj = Heaviside(model)
            obj.model = model;
            obj.dims = model.dims;
        end
        
        function generate_variables(obj, opts)
            switch class(obj.model)
              case "nosnoc.model.Pss"
                obj.generate_variables_from_pss(opts);
              case "nosnoc.model.Heaviside"
                obj.generate_variables_from_heaviside(opts);
              otherwise
                error('A heaviside DCS can only be created from a PSS or Heaviside step model.')
            end
        end

        function generate_equations(obj, opts)
            switch class(obj.model)
              case "nosnoc.model.Pss"
                obj.generate_equations_from_pss(opts);
              case "nosnoc.model.Heaviside"
                obj.generate_equations_from_heaviside(opts);
              otherwise
                error('A heaviside DCS can only be created from a PSS or Heaviside step model.')
            end
        end
    end

    methods(Access=private)
        function generate_variables_from_pss(obj, opts)
            import casadi.*
            dims = obj.dims;
            model = obj.model;

            % Check pss popluated S and c;
            if isempty(model.S)
                error("The PSS model must containt the switching matrix S and switching function c.");
            end
                
            % dimensions
            dims.n_alpha = sum(obj.dims.n_c_sys); % number of modes
            dims.n_lambda = dims.n_alpha;
            idx = 1;
            for ii = 1:dims.n_sys
                sys_idx{ii} = idx:(idx + obj.dims.n_c_sys(ii)-1);
                idx = idx + obj.dims.n_c_sys(ii);
                ii_str = num2str(ii);
                % define alpha (selection of a set valued step function)
                obj.alpha_sys{ii} = define_casadi_symbolic(opts.casadi_symbolic_mode,['alpha_' ii_str],obj.dims.n_c_sys(ii));
                obj.alpha = [obj.alpha;obj.alpha_sys{ii}];
                % define lambda_n_i (Lagrange multipler of alpha >= 0;)
                obj.lambda_n_sys{ii} = define_casadi_symbolic(opts.casadi_symbolic_mode,['lambda_n_' ii_str],obj.dims.n_c_sys(ii));
                obj.lambda_n = [obj.lambda_n;obj.lambda_n_sys{ii}];
                % define lambda_p_i (Lagrange multipler of alpha <= 1;)
                obj.lambda_p_sys{ii} = define_casadi_symbolic(opts.casadi_symbolic_mode,['lambda_p_' ii_str],obj.dims.n_c_sys(ii));
                obj.lambda_p = [obj.lambda_p;obj.lambda_p_sys{ii}];
            end

            % symbolic variables z = [alpha;obj.lambda_n;obj.lambda_p];
            obj.z_all = [obj.alpha;obj.lambda_n;obj.lambda_p;obj.model.z];
            obj.dims = dims;
        end

        function generate_variables_from_heaviside(obj, opts)
            import casadi.*
            dims = obj.dims;
            % dimensions

        end

        function generate_equations_from_pss(obj, opts)
            import casadi.*
            model = obj.model;
            dims = obj.dims;

            % generate theta_expressions
            for ii = 1:dims.n_sys
                theta_i = [];
                S_temp = model.S{ii};
                if opts.pss_lift_step_functions
                    % TODO implement automatic lifting
                else
                    for jj = 1:size(S_temp,1)
                        theta_jk = 1;
                        for kk = 1:size(S_temp,2)
                            % create multiaffine term
                            if S_temp(jj,kk) ~=0
                                theta_jk = theta_jk * (0.5*(1-S_temp(jj,kk))+S_temp(jj,kk)*obj.alpha_sys{ii}(kk));
                            end
                        end
                        theta_i = [theta_i;theta_jk];
                    end                    
                end
                obj.theta_expr_sys{ii} = theta_i;
            end

            obj.f_x = zeros(dims.n_x,1);
            for ii = 1:dims.n_sys
                obj.f_x = obj.f_x + model.F{ii}*obj.theta_expr_sys{ii};
            end

            g_switching = []; % collects switching function algebraic equations, 0 = c(x)-lambda_p+lambda_n
            lambda00_expr =[];
            for ii = 1:dims.n_sys
                % c_i(x) - (lambda_p_i-lambda_n_i)  = 0; for all i = 1,..., n_sys
                % lambda_n_i'*alpha_i  = 0; for all i = 1,..., n_sys
                % lambda_p_i'*(e-alpha_i)  = 0; for all i = 1,..., n_sys
                % lambda_n_i >= 0;    for all i = 1,..., n_sys
                % lambda_p_i >= 0;    for all i = 1,..., n_sys
                % alpha_i >= 0;     for all i = 1,..., n_sys
                g_switching = [g_switching;model.c{ii}-obj.lambda_p_sys{ii}+obj.lambda_n_sys{ii}];
                lambda00_expr = [lambda00_expr; -min(model.c{ii}, 0); max(model.c{ii},0)];
            end
            g_alg = [g_switching];

            obj.f_x_fun = Function('f_x', {model.x, model.z, obj.alpha, obj.lambda_n, obj.lambda_p, model.u, model.v_global, model.p}, {obj.f_x, model.f_q});
            obj.f_q_fun = Function('f_q', {model.x, model.z, obj.alpha, obj.lambda_n, obj.lambda_p, model.u, model.v_global, model.p}, {model.f_q});
            obj.g_z_fun = Function('g_z', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_z});
            obj.g_alg_fun = Function('g_alg', {model.x, model.z, obj.alpha, obj.lambda_n, obj.lambda_p, model.u, model.v_global, model.p}, {g_alg});
            obj.g_switching_fun = Function('g_switching', {model.x, model.z, obj.lambda_n, obj.lambda_p, model.v_global, model.p}, {g_switching});
            obj.lambda00_fun = Function('lambda00', {model.x, model.z, model.v_global, model.p_global}, {lambda00_expr});
            obj.g_path_fun = Function('g_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_path}); % TODO(@anton) do dependence checking for spliting the path constriants
            obj.g_comp_path_fun  = Function('g_comp_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_comp_path});
            obj.g_terminal_fun  = Function('g_terminal', {model.x, model.z, model.v_global, model.p_global}, {model.g_terminal});
            obj.f_q_T_fun = Function('f_q_T', {model.x, model.z, model.v_global, model.p}, {model.f_q_T});
            obj.f_lsq_x_fun = Function('f_lsq_x_fun',{model.x,model.x_ref,model.p},{model.f_lsq_x});
            obj.f_lsq_u_fun = Function('f_lsq_u_fun',{model.u,model.u_ref,model.p},{model.f_lsq_u});
            obj.f_lsq_T_fun = Function('f_lsq_T_fun',{model.x,model.x_ref_end,model.p_global},{model.f_lsq_T});
        end

        function generate_equations_from_heaviside(obj, opts)
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
