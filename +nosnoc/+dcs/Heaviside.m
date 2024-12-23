classdef Heaviside < nosnoc.dcs.Base
    properties
        alpha % casadi.SX|casadi.MX: Heaviside step variables
        alpha_sys % cell(casadi.SX|casadi.MX): Cell array with each cell containing alpha for linearly independent subsystem.
        lambda_n % casadi.SX|casadi.MX: Variables representing the negative part of the switching function.
        lambda_n_sys % cell(casadi.SX|casadi.MX): Cell array with each cell containing lambda_n for linearly independent subsystem.
        lambda_p % casadi.SX|casadi.MX: Variables representing the positive part of the switching function.
        lambda_p_sys % cell(casadi.SX|casadi.MX): Cell array with each cell containing lambda_p for linearly independent subsystem.
        theta_expr_sys % cell(casadi.SX|casadi.MX): Cell array with each cell containing the expression in terms of alpha of the convex multiplier theta

        z_all % casadi.SX|casadi.MX: Vector of all algorithmic algebraic variables.

        f_x % casadi.SX|casadi.MX: Expression for differential part of the DCS r.h.s: $\theta^\top F(x,u)%.
        g_alg % casadi.SX|casadi.MX: Expression for the algebraic part of the DCS. 
        
        dims % struct: Contains relevant dimensions.

        g_lp_stationarity_fun % casadi.Function: Function for the Heaviside LP stationarity for all switching functions.
        lambda00_fun % casadi.Function: Function mapping initial lambda_n and lambda_p from initial x.
    end

    properties(SetAccess=private)
        z_depends_on_alpha(1,1) logical = false; % boolean: defines whether the user algebraics z depend on alpha, as this requires some special handling. 
    end


    methods
        function obj = Heaviside(model)
            obj.model = model;
            obj.dims = model.dims;
        end
        
        function generate_variables(obj, opts)
            import casadi.*
            dims = obj.dims;
            model = obj.model;

            if ~ismember(class(obj.model), {'nosnoc.model.Pss', 'nosnoc.model.Heaviside'})
                nosnoc.error('wrong_model_type', 'A heaviside DCS can only be created from a PSS or Heaviside step model.')
            end

            switch class(obj.model)
              case "nosnoc.model.Pss"
                % Check pss popluated S and c;
                if isempty(model.S)
                    nosnoc.error('missing_S', "The PSS model must containt the lp_stationarity matrix S and lp_stationarity function c.");
                end

                dims.n_alpha = sum(obj.dims.n_c_sys); % number of modes
              case "nosnoc.model.Heaviside"
                obj.alpha = model.alpha;
                if all(size(model.g_z) ~= 0) && model.g_z.depends_on(obj.alpha)
                    obj.z_depends_on_alpha = true;
                end
            end
            dims.n_lambda = dims.n_alpha;

            idx = 1;
            for ii = 1:dims.n_sys
                sys_idx{ii} = idx:(idx + obj.dims.n_c_sys(ii)-1);
                idx = idx + obj.dims.n_c_sys(ii);
                ii_str = num2str(ii);
                % define alpha (selection of a set valued step function)
                switch class(obj.model)
                  case "nosnoc.model.Pss"
                    obj.alpha_sys{ii} = define_casadi_symbolic(opts.casadi_symbolic_mode,['alpha_' ii_str],obj.dims.n_c_sys(ii));
                    obj.alpha = [obj.alpha;obj.alpha_sys{ii}];
                  case "nosnoc.model.Heaviside"
                    obj.alpha_sys{ii} = obj.alpha(sys_idx{ii});
                end
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

        function generate_equations(obj, opts)
            import casadi.*
            model = obj.model;
            dims = obj.dims;

            switch class(obj.model)
              case "nosnoc.model.Pss"
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
              case "nosnoc.model.Heaviside"
                obj.f_x = model.f_x;
            end

            g_lp_stationarity = []; % collects lp_stationarity function algebraic equations, 0 = c(x)-lambda_p+lambda_n
            lambda00_expr =[];
            for ii = 1:dims.n_sys
                % c_i(x) - (lambda_p_i-lambda_n_i)  = 0; for all i = 1,..., n_sys
                % lambda_n_i'*alpha_i  = 0; for all i = 1,..., n_sys
                % lambda_p_i'*(e-alpha_i)  = 0; for all i = 1,..., n_sys
                % lambda_n_i >= 0;    for all i = 1,..., n_sys
                % lambda_p_i >= 0;    for all i = 1,..., n_sys
                % alpha_i >= 0;     for all i = 1,..., n_sys
                g_lp_stationarity = [g_lp_stationarity;model.c{ii}-obj.lambda_p_sys{ii}+obj.lambda_n_sys{ii}];
                lambda00_expr = [lambda00_expr; -min(model.c{ii}, 0); max(model.c{ii},0)];
            end
            g_alg = [g_lp_stationarity];
            obj.g_alg = g_alg;

            obj.f_x_fun = Function('f_x', {model.x, model.z, obj.alpha, obj.lambda_n, obj.lambda_p, model.u, model.v_global, model.p}, {obj.f_x, model.f_q});
            obj.f_q_fun = Function('f_q', {model.x, model.z, obj.alpha, obj.lambda_n, obj.lambda_p, model.u, model.v_global, model.p}, {model.f_q});
            obj.g_z_fun = Function('g_z', {model.x, obj.alpha, model.z, model.u, model.v_global, model.p}, {model.g_z}); % alpha taken from dcs rather than model as alpha only exists in Heaviside and not PSS models.
            obj.g_alg_fun = Function('g_alg', {model.x, model.z, obj.alpha, obj.lambda_n, obj.lambda_p, model.u, model.v_global, model.p}, {g_alg});
            obj.g_lp_stationarity_fun = Function('g_lp_stationarity', {model.x, model.z, obj.lambda_n, obj.lambda_p, model.v_global, model.p}, {g_lp_stationarity});
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
            var_list.lambda_p = obj.lambda_p;
            var_list.lambda_n = obj.lambda_n;
            var_list.alpha = obj.alpha;
            propgrp(2) = matlab.mixin.util.PropertyGroup(var_list, group_title);
        end
    end
end
