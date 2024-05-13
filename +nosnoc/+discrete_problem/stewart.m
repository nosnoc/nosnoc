classdef stewart < vdx.problems.Mpcc
    properties
        model
        dcs
        opts
    end

    methods
        function obj = stewart(dcs, opts)
            obj = obj@vdx.problems.Mpcc();
            obj.model = dcs.model
            obj.dcs = dcs;
            obj.opts = opts;
        end

        function create_variables(obj)
            dims = obj.dcs.dims;
            dcs = obj.dcs;
            model = obj.model;
            
            obj.p.rho_h_p = {{'rho_h_p',1}, 1};
            obj.p.T = {{'T',1}, opts.T};
            obj.p.p_global = {model.p_global, model.p_global_val};

            % other derived values
            h0 = opts.h;

            % 0d vars
            obj.w.v_global = {{'v_global',dims.n_v_global}, model.lbv_global, model.ubv_global, model.v0_global};

            % 1d vars
            obj.w.u(1:opts.N_stages) = {{'u', dims.n_u}, model.lbu, model.ubu, model.u0};
            obj.p.p_time_var(1:opts.N_stages) = {{'p_time_var', dims.n_p_time_var}, model.p_time_var_val};
            if obj.opts.use_speed_of_time_variables
                obj.w.sot(1:opts.N_stages) = {{'sot', 1}, opts.s_sot_min, opts.s_sot_max, opts.s_sot0};
            end

            % 2d vars
            if obj.opts.use_fesd
                obj.w.h(1:opts.N_stages,1:opts.N_finite_elements(1)) = {{'h', 1}, (1-opts.gamma_h)*h0, (1+opts.gamma_h)*h0, h0};
            end
            if (strcmp(obj.opts.step_equilibration,'linear')||..
                strcmp(obj.opts.step_equilibration,'linear_tanh')||...
                strcmp(obj.opts.step_equilibration,'linear_relaxed'))
                obj.w.B_max(1:opts.N_stages,2:opts.N_finite_elements) = {{'B_max', dims.n_c},-inf,inf};
                obj.w.pi_lambda(1:opts.N_stages,2:opts.N_finite_elements(1)) = {{'pi_lambda', dims.n_c},-inf,inf};
                obj.w.pi_c(1:opts.N_stages,2:opts.N_finite_elements(1)) = {{'pi_c', dims.n_c},-inf,inf};
                obj.w.lambda_lambda(1:opts.N_stages,2:opts.N_finite_elements(1)) = {{'lambda_lambda', dims.n_c},0,inf};
                obj.w.lambda_c(1:opts.N_stages,2:opts.N_finite_elements(1)) = {{'lambda_c', dims.n_c},0,inf};
                obj.w.eta(1:opts.N_stages,2:opts.N_finite_elements(1)) = {{'eta', dims.n_c},0,inf};
                obj.w.nu(1:opts.N_stages,2:opts.N_finite_elements(1)) = {{'nu', 1},0,inf};
            end

            % 3d vars
            obj.w.x(0,0,opts.n_s) = {{['x_0'], dims.n_x}};
            obj.w.x(1:opts.N_stages,1:opts.N_finite_elements(1),1:opts.n_s) = {{'x', dims.n_x}, model.lbx, model.ubx, model.x0};
            obj.w.z(1:opts.N_stages,1:opts.N_finite_elements(1),1:opts.n_s) = {{'z', dims.n_z}, model.lbz, model.ubz, model.z0};
            obj.w.lambda(0,0,opts.n_s) = {{['lambda_0'], dims.n_lambda},0,inf};
            obj.w.lambda(1:opts.N_stages,1:opts.N_finite_elements(1),1:opts.n_s) = {{'lambda', dims.n_lambda},0, inf};
            obj.w.theta(1:opts.N_stages,1:opts.N_finite_elements(1),1:opts.n_s) = {{'theta', dims.n_theta},0, 1};
            obj.w.mu(1:opts.N_stages,1:opts.N_finite_elements(1),1:opts.n_s) = {{'mu', dims.n_mu},0,inf};
        end

        function forward_sim_constraints(obj)
            import casadi.*
            model = obj.model;
            opts = obj.opts;
            if obj.opts.use_fesd
                t_stage = obj.p.T()/opts.N_stages;
                h0 = obj.p.T().val/(opts.N_stages*opts.N_finite_elements(1));
            else
                h0 = obj.p.T().val/(opts.N_stages*opts.N_finite_elements(1));
            end
            
            v_global = obj.w.v_global();
            p_global = obj.p.p_global();

            x_0 = obj.w.x(0,0,opts.n_s);
            lambda_0 = obj.w.lambda(0,0,opts.n_s);
            
            x_prev = obj.w.x(0,0,opts.n_s);
            for ii=1:opts.N_stages
                ui = obj.w.u(ii);
                p_stage = obj.p.p_time_var(ii);
                p = [p_global;p_stage];
                if obj.opts.use_speed_of_time_variables
                    s_sot = obj.w.sot(ii);
                else
                    s_sot = 1;
                end
                
                sum_h = 0;
                for jj=1:opts.N_finite_elements(ii)
                    if obj.opts.use_fesd
                        h = obj.w.h(ii,jj);
                        sum_h = sum_h + h;
                    else
                        h = h0;
                    end
                    for kk=1:opts.n_s
                        x_ijk = obj.w.x(ii,jj,kk);
                        z_ijk = obj.w.z(ii,jj,kk);
                        lambda_ijk = obj.w.lambda(ii,jj,kk);
                        theta_ijk = obj.w.theta(ii,jj,kk);
                        mu_ijk = obj.w.mu(ii,jj,kk);
                        
                        
                        fj = s_sot*dcs.f_x_fun(x_ijk, z_ijk, lambda_ijk, theta_ijk, mu_ijk, ui, v_global, p);
                        qj = s_sot*dcs.f_q_fun(x_ijk, z_ijk, lambda_ijk, theta_ijk, mu_ijk, ui, v_global, p);
                        xk = opts.C_irk(1, kk+1) * x_prev;
                        for rr=1:opts.n_s % TODO(@anton) handle other modes.
                            x_ijr = obj.w.x(ii,jj,rr);
                            xk = xk + opts.C_irk(rr+1, kk+1) * x_ijr;
                        end
                        obj.g.dynamics(ii,jj,kk) = {h * fj - xk};
                        obj.g.z(ii,jj,kk) = {dcs.g_z_fun(x_ijk, z_ijk, ui, v_global, p)};
                        obj.g.path(ii,jj,kk) = {dcs.g_path_fun(x_ijk, z_ijk, ui, v_global, p_global, p_stage), model.lbg_path, model.ubg_path};
                        obj.g.algebraic(ii,jj,kk) = {dcs.g_alg_fun(x_ijk, z_ijk, lambda_ijk, theta_ijk, mu_ijk, ui, v_global, p)};
                        % also integrate the objective
                        obj.f = obj.f + opts.B_irk(kk+1)*h*qj;
                    end
                    x_prev = obj.w.x(ii,jj,opts.n_s);
                end
                if obj.opts.use_fesd
                    obj.g.sum_h(ii) = {t_stage-sum_h};
                end
                if obj.opts.use_speed_of_time_variables
                    x_end = obj.w.x(ii,opts.N_finite_elements(end),opts.n_s);
                    x_start = obj.w.x(0,0,opts.n_s);
                    obj.g.g_equidistant_grid(ii) = {(x_end(end)-x_start(end)) - t_stage*ii};
                end
            end

            % Terminal cost
            obj.f = obj.f + model.f_q_T_fun(obj.w.x(ii,jj,kk), obj.w.z(ii,jj,kk), v_global, p_global);

            % Terminal constraint
            obj.g.terminal = {model.g_terminal_fun(obj.w.x(ii,jj,kk), obj.w.z(ii,jj,kk), v_global, p_global), model.lbg_terminal, model.ubg_terminal};
        end

        function generate_complementarities(obj)
            import casadi.*
            opts = obj.opts;
            dcs = obj.dcs;
            model = obj.model;
            % Do Cross-Complementarity

            if opts.use_fesd
                switch opts.cross_comp_mode
                  case CrossCompMode.STAGE_STAGE
                    for ii=1:opts.N_stages
                        for jj=1:opts.N_finite_elements(ii);
                            Gij = {};
                            Hij = {};
                            for kk=1:opts.n_s
                                lambda_ijk = obj.w.theta(ii,jj,kk);
                                for rr=1:opts.n_s
                                    theta_ijr = obj.w.theta(ii,jj,rr);
                                    
                                    Gij = vertcat(Gij, lambda_ijk);
                                    Hij = vertcat(Hij, theta_ijr);
                                end
                            end
                            obj.G.cross_comp(ii,jj) = {vertcat(Gij{:})};
                            obj.H.cross_comp(ii,jj) = {vertcat(Hij{:})};
                            lambda_prev = obj.w.lambda(ii,jj,opts.n_s);
                        end
                    end
                  case CrossCompMode.FE_STAGE
                    lambda_prev = obj.w.lambda(0,0,opts.n_s);
                    % Do cross comp for distance with lambda
                    for ii=1:opts.N_stages
                        for jj=1:opts.N_finite_elements(ii);
                            sum_lambda = lambda_prev + sum1(obj.w.lambda(ii,jj,:));
                            Gij = {};
                            Hij = {};
                            for kk=1:opts.n_s
                                theta_ijk = obj.w.theta(ii,jj,kk);
                                
                                Gij = vertcat(Gij, sum_lambda);
                                Hij = vertcat(Hij, theta_ijk);
                            end
                            obj.G.cross_comp(ii,jj) = {vertcat(Gij{:})};
                            obj.H.cross_comp(ii,jj) = {vertcat(Hij{:})};
                            lambda_prev = obj.w.lambda(ii,jj,opts.n_s);
                        end
                    end
                  case CrossCompMode.STAGE_FE
                    lambda_prev = obj.w.lambda(0,0,opts.n_s);
                    % Do cross comp for distance with lambda
                    for ii=1:opts.N_stages
                        for jj=1:opts.N_finite_elements(ii);
                            sum_theta = sum1(obj.w.theta(ii,jj,:));
                            Gij = {lambda_prev};
                            Hij = {sum_theta};
                            for kk=1:opts.n_s
                                lambda_ijk = obj.w.lambda(ii,jj,kk);
                                
                                Gij = vertcat(Gij, lambda_ijk);
                                Hij = vertcat(Hij, sum_theta);
                            end
                            obj.G.cross_comp(ii,jj) = {vertcat(Gij{:})};
                            obj.H.cross_comp(ii,jj) = {vertcat(Hij{:})};
                            lambda_prev = obj.w.lambda(ii,jj,opts.n_s);
                        end
                    end
                  case CrossCompMode.FE_FE
                    lambda_prev = obj.w.lambda(0,0,opts.n_s);
                    % Do cross comp for distance with lambda
                    for ii=1:opts.N_stages
                        for jj=1:opts.N_finite_elements(ii);
                            Gij = [];
                            Hij = [];
                            for kk=1:opts.n_s
                                lambda_ijk = obj.w.lambda(ii,jj,kk);
                                theta_ijk = obj.w.theta(ii,jj,kk);
                                
                                Gij = Gij + lambda_ijk;
                                Hij = Hij + theta_ijk;
                            end
                            obj.G.cross_comp(ii,jj) = {Gij};
                            obj.H.cross_comp(ii,jj) = {Hij};
                            lambda_prev = obj.w.lambda(ii,jj,opts.n_s);
                        end
                    end
                end
            else
                for ii=1:opts.N_stages
                    for jj=1:opts.N_finite_elements(ii);
                        Gij = {};
                        Hij = {};
                        for kk=1:opts.n_s
                            lambda_ijk = obj.w.theta(ii,jj,kk);
                            theta_ijk = obj.w.theta(ii,jj,kk);
                            
                            obj.G.standard_comp(ii,jj, kk) = {lambda_ijk};
                            obj.H.standard_comp(ii,jj, kk) = {theta_ijk};
                        end
                    end
                end
            end
        end

        function step_equilibration(obj)

        end

        function populate_problem(obj)
            obj.create_variables();
            obj.forward_sim_constraints();
            obj.generate_complementarities();
            obj.step_equilibration();

            obj.populated = true;
        end
        
        function create_solver(obj, solver_options, plugin)
            if ~obj.populated
                obj.populate_problem()
            end
            
            if ~exist('plugin')
                plugin = 'scholtes_ineq';
            end
            obj.w.sort_by_index();
            obj.g.sort_by_index();

            solver_options.assume_lower_bounds = true;
            
            create_solver@vdx.problems.Mpcc(obj, solver_options, plugin);
        end
    end
end
