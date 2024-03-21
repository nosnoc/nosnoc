classdef InclusionProblem < vdx.Problem
    properties (Access=public)
        data
        opts
        eta_vec
        eta_fun
        G_fun
        H_fun

        c_fun

        comp_res_fun
    end
    
    methods (Access=public)
        function obj = InclusionProblem(data, opts)
            obj@vdx.Problem(); % call parent constructor to prepare
            obj.data = data;
            obj.opts = opts;

            if ~isfield(obj.opts, 'use_fesd')
                obj.opts.use_fesd = true;
            end
            if ~isfield(obj.opts, 'step_eq')
                obj.opts.step_eq = 'heuristic_mean';
            end
            if ~obj.opts.use_fesd
                obj.opts.step_eq = 'none';
            end
            if ~isfield(obj.opts,'comp_scale')
                obj.opts.comp_scale = 1;
            end
            if ~isfield(obj.opts,'elastic_ell_inf')
                obj.opts.elastic_ell_inf = 0;
            end
            if ~isfield(obj.data, 'g_path')
                obj.data.g_path = [];
            end
            if ~isfield(obj.data, 'lbg_path')
                obj.data.lbg_path = zeros(size(obj.data.g_path));
            end
            if ~isfield(obj.data, 'ubg_path')
                obj.data.ubg_path = zeros(size(obj.data.g_path));
            end
            if ~isfield(obj.data, 'partial_proj_matrix')
                obj.data.partial_proj_matrix = eye(length(obj.data.x));
            end
            if ~isfield(obj.opts,'time_dependent')
                obj.opts.time_dependent = false;
            end
            % get dimensions
            n_u = length(data.u);
            n_c = length(data.c);
            n_x = length(data.x);

            if obj.opts.elastic_ell_inf
                obj.w.s_elastic(1) = {{'s_elastic',1},0,inf,0};
            end
            obj.p.sigma(1) = {{'sigma',1},0,inf,0};
            obj.p.gamma_h(1) = {{'gamma_h',1},0,inf,1e-1};
            obj.p.T(1) = {{'T',1},0,inf,data.T};

            % other derived values
            t_stage = data.T/data.N_stages;
            h0 = t_stage/data.N_fe;
            
            obj.w.x(0,0,data.n_s) = {{['x_0'], n_x}};
            obj.w.lambda(0,0,data.n_s) = {{['lambda_0'], n_c},0,0};
            for ii=1:data.N_stages
                obj.w.u(ii) = {{['u_' num2str(ii)], n_u}, data.lbu, data.ubu, data.u0};
                if obj.opts.time_dependent
                    obj.w.sot(ii) = {{['sot_' num2str(ii)], 1}, 0, 100, 1};
                end
                for jj=1:data.N_fe
                    if obj.opts.use_fesd
                        obj.w.h(ii,jj) = {{['h_' num2str(ii) '_' num2str(jj)], 1}, 0, 2*h0, h0};
                    end
                    if (strcmp(obj.opts.step_eq,'linear')||...
                        strcmp(obj.opts.step_eq,'linear_tanh')||...
                        strcmp(obj.opts.step_eq,'linear_relaxed')) && jj > 1
                        obj.w.B_max(ii,jj) ={{['B_max_' num2str(ii) '_' num2str(jj)], n_c},-inf,inf};
                        obj.w.pi_lambda(ii,jj) ={{['pi_lambda_' num2str(ii) '_' num2str(jj)], n_c},-inf,inf};
                        obj.w.pi_c(ii,jj) ={{['pi_c_' num2str(ii) '_' num2str(jj)], n_c},-inf,inf};
                        obj.w.lambda_lambda(ii,jj) ={{['lambda_lambda_' num2str(ii) '_' num2str(jj)], n_c},0,inf};
                        obj.w.lambda_c(ii,jj) ={{['lambda_c_' num2str(ii) '_' num2str(jj)], n_c},0,inf};
                        obj.w.eta(ii,jj) ={{['eta_' num2str(ii) '_' num2str(jj)], n_c},0,inf};
                        obj.w.nu(ii,jj) ={{['nu_' num2str(ii) '_' num2str(jj)], 1},0,inf};
                    end
                    for kk=1:data.n_s
                        obj.w.x(ii,jj,kk) = {{['x_' num2str(ii) '_' num2str(jj) '_' num2str(kk)], n_x}, data.lbx, data.ubx, data.x0};
                        obj.w.lambda(ii,jj,kk) = {{['lambda_' num2str(ii) '_' num2str(jj) '_' num2str(kk)], n_c},0,inf};
                    end
                end
            end
        end

        function generate_constraints(obj)
            import casadi.*
            [B, C, D, tau_root] = generate_butcher_tableu_integral(obj.data.n_s, obj.data.irk_scheme);
            %[~, B, C, ~] = generate_butcher_tableu(obj.data.n_s, 'RADAU_IIA')
            n_x = length(obj.data.x);
            n_u = length(obj.data.u);
            n_c = length(obj.data.c);

            % other derived values
            if obj.opts.use_fesd
                t_stage = obj.p.T(1)/obj.data.N_stages;
                h0 = obj.p.T(1).init/(obj.data.N_stages*obj.data.N_fe);
            else
                h0 = obj.p.T(1).init/(obj.data.N_stages*obj.data.N_fe);
            end
            
            if obj.opts.elastic_ell_inf
                sigma = obj.w.s_elastic(1);
            else
                sigma = obj.p.sigma(1);
            end
            % Define functions from obj.data
            lambda = SX.sym('lambda', n_c);
            
            nabla_c = obj.data.c.jacobian(obj.data.x)';

            E = obj.data.partial_proj_matrix;
            
            f_x = obj.data.f_x + E*nabla_c*lambda;
            f_x = f_x;

            g_path_fun = Function('g_path_fun', {obj.data.x}, {obj.data.g_path});
            f_x_fun = Function('f_x_fun', {obj.data.x,obj.data.u,lambda}, {f_x});
            f_q_fun = Function('q_fun', {obj.data.x,obj.data.u}, {obj.data.f_q});
            f_q_T_fun = Function('q_fun', {obj.data.x}, {obj.data.f_q_T});
            c_fun = Function('c_fun', {obj.data.x}, {obj.data.c});
            obj.c_fun = c_fun;
            c_dot_fun = Function('c_dot_fun', {obj.data.x,obj.data.u,lambda}, {nabla_c'*obj.data.f_x});
            x_prev = obj.w.x(0,0,obj.data.n_s);
            for ii=1:obj.data.N_stages
                ui = obj.w.u(ii);
                if obj.opts.time_dependent
                    s_sot = obj.w.sot(ii);
                else
                    s_sot = 1;
                end
                sum_h = 0;
                for jj=1:obj.data.N_fe
                    if obj.opts.use_fesd
                        h = obj.w.h(ii,jj);
                        sum_h = sum_h + h;
                    else
                        h = h0;
                    end
                    for kk=1:obj.data.n_s
                        x_ijk = obj.w.x(ii,jj,kk);
                        lambda_ijk = obj.w.lambda(ii,jj,kk);
                        fj = s_sot*f_x_fun(x_ijk,ui,lambda_ijk);
                        qj = s_sot*f_q_fun(x_ijk,ui);
                        xk = C(1, kk+1) * x_prev;
                        for rr=1:obj.data.n_s
                            x_ijr = obj.w.x(ii,jj,rr);
                            xk = xk + C(rr+1, kk+1) * x_ijr;
                        end
                        obj.g.dynamics(ii,jj,kk) = {h * fj - xk};
                        % also add non-negativity constraint on c
                        obj.g.c_nonnegative(ii,jj,kk) = {c_fun(x_ijk), 0, inf};
                        %add path constraints
                        obj.g.g_path(ii,jj,kk) = {g_path_fun(x_ijk), obj.data.lbg_path, obj.data.ubg_path};
                        % also integrate the objective
                        obj.f = obj.f + B(kk+1)*h*qj;
                    end
                    x_prev = obj.w.x(ii,jj,obj.data.n_s);
                end
                if obj.opts.use_fesd
                    obj.g.sum_h(ii) = {t_stage-sum_h};
                end
                if obj.opts.time_dependent
                    x_end = obj.w.x(ii,obj.data.N_fe,obj.data.n_s);
                    x_start = obj.w.x(0,0,obj.data.n_s);
                    obj.g.g_equidistant_grid(ii) = {(x_end(end)-x_start(end)) - t_stage*ii};
                end
            end

            % Terminal cost
            obj.f = obj.f + f_q_T_fun(obj.w.x(ii,jj,kk));


            % Terminal constraint
            if isfield(obj.data, 'g_T')
                g_T_fun = Function('g_T_fun', {obj.data.x}, {obj.data.g_T});
                obj.g.terminal(0) = {g_T_fun(obj.w.x(ii,jj,kk))}; % TODO(@anton) assume equality for now
            end
            
            % Do Cross-Complementarity
            x_prev = obj.w.x(0,0,obj.data.n_s);
            lambda_prev = obj.w.lambda(0,0,obj.data.n_s);
            G = [];
            H = [];
            for ii=1:obj.data.N_stages
                for jj=1:obj.data.N_fe
                    if obj.opts.use_fesd
                        Gij = c_fun(x_prev);
                        Hij = lambda_prev;
                        for kk=1:obj.data.n_s
                            x_ijk = obj.w.x(ii,jj,kk);
                            lambda_ijk = obj.w.lambda(ii,jj,kk);
                            Gij = Gij + c_fun(x_ijk);
                            Hij = Hij + lambda_ijk;
                        end

                        G = [G;Gij];
                        H = [H;Hij];
                        
                        obj.g.complementarity(ii,jj) = {obj.opts.comp_scale*(Gij.*Hij - sigma), -inf, 0};
                        x_prev = obj.w.x(ii,jj,obj.data.n_s);
                        lambda_prev = obj.w.lambda(ii,jj,obj.data.n_s);
                    else
                        Gij = [];
                        Hij = [];
                        for kk=1:obj.data.n_s
                            x_ijk = obj.w.x(ii,jj,kk);
                            lambda_ijk = obj.w.lambda(ii,jj,kk);
                            Gij = [Gij;c_fun(x_ijk)];
                            Hij = [Hij;lambda_ijk];
                        end
                        G = [G;Gij];
                        H = [H;Hij];
                        obj.g.complementarity(ii,jj) = {obj.opts.comp_scale*(Gij.*Hij - sigma), -inf, 0};
                    end
                end
            end
            obj.G_fun = Function('G_fun',{obj.w.w},{G});
            obj.H_fun = Function('H_fun',{obj.w.w},{H});

            
            eta_vec = [];
            for ii=1:obj.data.N_stages
                for jj=2:obj.data.N_fe
                    sigma_c_B = 0;
                    sigma_lam_B = 0;
                    for kk=1:obj.data.n_s
                        sigma_c_B = sigma_c_B + c_fun(obj.w.x(ii,jj-1,kk));
                        sigma_lam_B = sigma_lam_B + obj.w.lambda(ii,jj-1,kk);
                    end
                    sigma_c_F = 0;
                    sigma_lam_F = 0;
                    for kk=1:obj.data.n_s
                        sigma_c_F = sigma_c_F + c_fun(obj.w.x(ii,jj,kk));
                        sigma_lam_F = sigma_lam_F + obj.w.lambda(ii,jj,kk);
                    end

                    pi_c = sigma_c_B .* sigma_c_F;
                    pi_lam = sigma_lam_B .* sigma_lam_F;
                    nu = pi_c + pi_lam;
                    eta = 1;
                    for jjj=1:length(nu)
                        eta = eta*nu(jjj);
                    end
                    eta_vec = [eta_vec;eta];
                end
            end
            obj.eta_fun = Function('eta_fun', {obj.w.w}, {eta_vec});
            switch obj.opts.step_eq
              case 'heuristic_mean'
                for ii=1:obj.data.N_stages
                    for jj=1:obj.data.N_fe
                        obj.f = obj.f + obj.p.gamma_h(1)*(h0-obj.w.h(ii,jj))^2;
                    end
                end
              case 'direct'
                eta_vec = [];
                for ii=1:obj.data.N_stages
                    for jj=2:obj.data.N_fe
                        sigma_c_B = 0;
                        sigma_lam_B = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_B = sigma_c_B + c_fun(obj.w.x(ii,jj-1,kk));
                            sigma_lam_B = sigma_lam_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lam_F = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_F = sigma_c_F + c_fun(obj.w.x(ii,jj,kk));
                            sigma_lam_F = sigma_lam_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lam_B .* sigma_lam_F;
                        nu = pi_c + pi_lam;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        obj.g.step_equilibration(ii,jj) = {eta*delta_h, 0, 0};
                    end
                end
                obj.eta_fun = Function('eta_fun', {obj.w.w}, {eta_vec});
              case 'direct_homotopy'
                eta_vec = [];
                for ii=1:obj.data.N_stages
                    for jj=2:obj.data.N_fe
                        sigma_c_B = 0;
                        sigma_lam_B = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_B = sigma_c_B + c_fun(obj.w.x(ii,jj-1,kk));
                            sigma_lam_B = sigma_lam_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lam_F = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_F = sigma_c_F + c_fun(obj.w.x(ii,jj,kk));
                            sigma_lam_F = sigma_lam_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lam_B .* sigma_lam_F;
                        nu = pi_c + pi_lam;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        obj.eta_vec = eta_vec;
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        homotopy_eq = [eta*delta_h - sigma;eta*delta_h + sigma];
                        obj.g.step_equilibration(ii,jj) = {homotopy_eq, [-inf;0], [0;inf]};
                    end
                end
                obj.eta_fun = Function('eta_fun', {obj.w.w}, {eta_vec});
              case 'direct_homotopy_with_penalty'
                eta_vec = [];
                for ii=1:obj.data.N_stages
                    for jj=2:obj.data.N_fe
                        sigma_c_B = 0;
                        sigma_lam_B = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_B = sigma_c_B + c_fun(obj.w.x(ii,jj-1,kk));
                            sigma_lam_B = sigma_lam_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lam_F = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_F = sigma_c_F + c_fun(obj.w.x(ii,jj,kk));
                            sigma_lam_F = sigma_lam_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lam_B .* sigma_lam_F;
                        nu = pi_c + pi_lam;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        obj.eta_vec = eta_vec;
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        homotopy_eq = [eta*delta_h - sigma;eta*delta_h + sigma];
                        obj.g.step_equilibration(ii,jj) = {homotopy_eq, [-inf;0], [0;inf]};
                    end
                end
                obj.eta_fun = Function('eta_fun', {obj.w.w}, {eta_vec});
                for ii=1:obj.data.N_stages
                    for jj=1:obj.data.N_fe
                        obj.f = obj.f + obj.p.gamma_h(1)*(h0-obj.w.h(ii,jj))^2;
                    end
                end
              case 'direct_fix_pathological'
                %TODO(@anton)
              case 'linear_tanh'
                for ii=1:obj.data.N_stages
                    for jj=2:obj.data.N_fe
                        sigma_c_b = 0;
                        sigma_lambda_b = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_b = sigma_c_b + c_fun(obj.w.x(ii,jj-1,kk));
                            sigma_lambda_b = sigma_lambda_b + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_f = 0;
                        sigma_lambda_f = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_f = sigma_c_f + c_fun(obj.w.x(ii,jj,kk));
                            sigma_lambda_f = sigma_lambda_f + obj.w.lambda(ii,jj,kk);
                        end

                        
                        % todo ideally we output G and H instead of doing all of the stuff here but ok.
                        lambda_c = obj.w.lambda_c(ii,jj);
                        lambda_lambda = obj.w.lambda_lambda(ii,jj);
                        B_max = obj.w.B_max(ii,jj);
                        pi_c = obj.w.pi_c(ii,jj);
                        pi_lambda = obj.w.pi_lambda(ii,jj);
                        eta = obj.w.eta(ii,jj);
                        nu = obj.w.nu(ii,jj);

                        obj.g.pi_c_or(ii,jj) = {[pi_c-sigma_c_f;pi_c-sigma_c_b;sigma_c_f+sigma_c_b-pi_c],0,inf};
                        obj.g.pi_lambda_or(ii,jj) = {[pi_lambda-sigma_lambda_f;pi_lambda-sigma_lambda_b;sigma_lambda_f+sigma_lambda_b-pi_lambda],0,inf};

                        % kkt conditions for min B, B>=sigmaB, B>=sigmaF
                        kkt_max = [1-lambda_lambda-lambda_c;
                            B_max-pi_c;
                            B_max-pi_lambda;
                            (B_max-pi_c).*lambda_c - sigma;
                            (B_max-pi_lambda).*lambda_lambda - sigma];
                        obj.g.kkt_max(ii,jj) = {kkt_max,
                            [0*ones(n_c,1);0*ones(n_c,1);0*ones(n_c,1);-inf*ones(n_c,1);-inf*ones(n_c,1)],
                            [0*ones(n_c,1);inf*ones(n_c,1);inf*ones(n_c,1);0*ones(n_c,1);0*ones(n_c,1)]};

                        % eta calculation
                        eta_const = [eta-pi_lambda;eta-pi_c;eta-pi_lambda-pi_c+B_max];
                        obj.g.eta_const(ii,jj) = {eta_const,
                            [-inf*ones(n_c,1);-inf*ones(n_c,1);zeros(n_c,1)],
                            [zeros(n_c,1);zeros(n_c,1);inf*ones(n_c,1)]};

                        obj.g.nu_or(ii,jj) = {[nu-eta;sum(eta)-nu],0,inf};

                        % the actual step eq conditions
                        %M = 1e5;
                        M=100*t_stage;
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        step_equilibration = [delta_h + tanh(1000*nu)*M;
                            delta_h - tanh(1000*nu)*M];
                        obj.g.step_equilibration(ii,jj) = {step_equilibration,[0;-inf],[inf;0]};
                        
                    end
                end
              case 'linear'
                for ii=1:obj.data.N_stages
                    for jj=2:obj.data.N_fe
                        sigma_c_b = 0;
                        sigma_lambda_b = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_b = sigma_c_b + c_fun(obj.w.x(ii,jj-1,kk));
                            sigma_lambda_b = sigma_lambda_b + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_f = 0;
                        sigma_lambda_f = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_f = sigma_c_f + c_fun(obj.w.x(ii,jj,kk));
                            sigma_lambda_f = sigma_lambda_f + obj.w.lambda(ii,jj,kk);
                        end

                        
                        % todo ideally we output G and H instead of doing all of the stuff here but ok.
                        lambda_c = obj.w.lambda_c(ii,jj);
                        lambda_lambda = obj.w.lambda_lambda(ii,jj);
                        B_max = obj.w.B_max(ii,jj);
                        pi_c = obj.w.pi_c(ii,jj);
                        pi_lambda = obj.w.pi_lambda(ii,jj);
                        eta = obj.w.eta(ii,jj);
                        nu = obj.w.nu(ii,jj);

                        obj.g.pi_c_or(ii,jj) = {[pi_c-sigma_c_f;pi_c-sigma_c_b;sigma_c_f+sigma_c_b-pi_c],0,inf};
                        obj.g.pi_lambda_or(ii,jj) = {[pi_lambda-sigma_lambda_f;pi_lambda-sigma_lambda_b;sigma_lambda_f+sigma_lambda_b-pi_lambda],0,inf};

                        % kkt conditions for min B, B>=sigmaB, B>=sigmaF
                        kkt_max = [1-lambda_lambda-lambda_c;
                            B_max-pi_c;
                            B_max-pi_lambda;
                            (B_max-pi_c).*lambda_c - sigma;
                            (B_max-pi_lambda).*lambda_lambda - sigma];
                        obj.g.kkt_max(ii,jj) = {kkt_max,
                            [0*ones(n_c,1);0*ones(n_c,1);0*ones(n_c,1);-inf*ones(n_c,1);-inf*ones(n_c,1)],
                            [0*ones(n_c,1);inf*ones(n_c,1);inf*ones(n_c,1);0*ones(n_c,1);0*ones(n_c,1)]};

                        % eta calculation
                        eta_const = [eta-pi_lambda;eta-pi_c;eta-pi_lambda-pi_c+B_max];
                        obj.g.eta_const(ii,jj) = {eta_const,
                            [-inf*ones(n_c,1);-inf*ones(n_c,1);zeros(n_c,1)],
                            [zeros(n_c,1);zeros(n_c,1);inf*ones(n_c,1)]};

                        obj.g.nu_or(ii,jj) = {[nu-eta;sum(eta)-nu],0,inf};

                        % the actual step eq conditions
                        %M = 1e5;
                        M=t_stage;
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        step_equilibration = [delta_h + (1/h0)*nu*M;
                            delta_h - (1/h0)*nu*M];
                        obj.g.step_equilibration(ii,jj) = {step_equilibration,[0;-inf],[inf;0]};
                        
                    end
                end
              case 'linear_relaxed'
                for ii=1:obj.data.N_stages
                    for jj=2:obj.data.N_fe
                        sigma_c_b = 0;
                        sigma_lambda_b = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_b = sigma_c_b + c_fun(obj.w.x(ii,jj-1,kk));
                            sigma_lambda_b = sigma_lambda_b + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_f = 0;
                        sigma_lambda_f = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_f = sigma_c_f + c_fun(obj.w.x(ii,jj,kk));
                            sigma_lambda_f = sigma_lambda_f + obj.w.lambda(ii,jj,kk);
                        end

                        
                        % todo ideally we output G and H instead of doing all of the stuff here but ok.
                        lambda_c = obj.w.lambda_c(ii,jj);
                        lambda_lambda = obj.w.lambda_lambda(ii,jj);
                        B_max = obj.w.B_max(ii,jj);
                        pi_c = obj.w.pi_c(ii,jj);
                        pi_lambda = obj.w.pi_lambda(ii,jj);
                        eta = obj.w.eta(ii,jj);
                        nu = obj.w.nu(ii,jj);

                        obj.g.pi_c_or(ii,jj) = {[pi_c-sigma_c_f;pi_c-sigma_c_b;sigma_c_f+sigma_c_b-pi_c],0,inf};
                        obj.g.pi_lambda_or(ii,jj) = {[pi_lambda-sigma_lambda_f;pi_lambda-sigma_lambda_b;sigma_lambda_f+sigma_lambda_b-pi_lambda],0,inf};

                        % kkt conditions for min B, B>=sigmaB, B>=sigmaF
                        kkt_max = [1-lambda_lambda-lambda_c;
                            B_max-pi_c;
                            B_max-pi_lambda;
                            (B_max-pi_c).*lambda_c - sigma;
                            (B_max-pi_lambda).*lambda_lambda - sigma];
                        obj.g.kkt_max(ii,jj) = {kkt_max,
                            [0*ones(n_c,1);0*ones(n_c,1);0*ones(n_c,1);-inf*ones(n_c,1);-inf*ones(n_c,1)],
                            [0*ones(n_c,1);inf*ones(n_c,1);inf*ones(n_c,1);0*ones(n_c,1);0*ones(n_c,1)]};

                        % eta calculation
                        eta_const = [eta-pi_lambda;eta-pi_c;eta-pi_lambda-pi_c+B_max];
                        obj.g.eta_const(ii,jj) = {eta_const,
                            [-inf*ones(n_c,1);-inf*ones(n_c,1);zeros(n_c,1)],
                            [zeros(n_c,1);zeros(n_c,1);inf*ones(n_c,1)]};

                        obj.g.nu_or(ii,jj) = {[nu-eta;sum(eta)-nu],0,inf};

                        % the actual step eq conditions
                        %M = 1e5;
                        obj.w.step_eq_slack(ii,jj) = {{['s_step_eq_' num2str(ii) '_' num2str(jj)],1},0,inf};
                        slack = obj.w.step_eq_slack(ii,jj);
                        M=100*t_stage;
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        step_equilibration = [delta_h + nu*M+slack;
                            delta_h - nu*M-slack];
                        obj.g.step_equilibration(ii,jj) = {step_equilibration,[0;-inf],[inf;0]};
                        obj.f = obj.f + (1/obj.p.sigma(1))*slack;
                    end
                end
              case 'none'
                % Do nothing
            end

            if obj.opts.elastic_ell_inf
                obj.f = obj.p.sigma(1)*obj.f + obj.w.s_elastic(1);
            end

            obj.comp_res_fun = Function('comp_res', {obj.w.w, obj.p.w}, {max(G.*H)});
        end

        % function stats = solve(obj)
        %     nabla_c = obj.data.c.jacobian(obj.data.x)';
        %     % Create initial cdot to populate the initial lambda parameter
        %     c_dot_fun = Function('c_dot_fun', {obj.data.x,obj.data.u,lambda}, {nabla_c'*obj.data.f_x});
        %     cdot0 = c_dot_fun(prob.w.x(0,0,data.n_s).init,)
        %     stats = solve@vdx.Problem(obj);
        % end
    end
end
