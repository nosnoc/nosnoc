classdef Gcs < vdx.problems.Mpcc
    properties
        model
        dcs
        opts

        populated = false
        sorted = false
    end

    methods
        function obj = Gcs(dcs, opts)
            obj = obj@vdx.problems.Mpcc();
            obj.model = dcs.model;
            obj.dcs = dcs;
            obj.opts = opts;
        end

        function create_variables(obj)
            dims = obj.dcs.dims;
            dcs = obj.dcs;
            model = obj.model;
            opts = obj.opts;

            % Parameters
            obj.p.rho_h_p = {{'rho_h_p',1}, 1};
            obj.p.rho_terminal_p = {{'rho_terminal_p',1}, 1};
            obj.p.T = {{'T',1}, opts.T};
            obj.p.p_global = {model.p_global, model.p_global_val};


            % 0d vars: Variables which only exist once globally. 
            obj.w.v_global = {{'v_global',dims.n_v_global}, model.lbv_global, model.ubv_global, model.v0_global};
            if opts.use_speed_of_time_variables && ~opts.local_speed_of_time_variable
                obj.w.sot = {{'sot', 1}, opts.s_sot_min, opts.s_sot_max, opts.s_sot0};
            end
            if opts.time_optimal_problem
                obj.w.T_final = {{'T_final', 1}, opts.T_final_min, opts.T_final_max, opts.T};
                obj.f = obj.f + obj.w.T_final();
            end

            % 1d vars: Variables that are defined per control stage
            obj.w.u(1:opts.N_stages) = {{'u', dims.n_u}, model.lbu, model.ubu, model.u0};
            obj.p.p_time_var(1:opts.N_stages) = {{'p_time_var', dims.n_p_time_var}, model.p_time_var_val};
            if opts.use_speed_of_time_variables && opts.local_speed_of_time_variable
                obj.w.sot(1:opts.N_stages) = {{'sot', 1}, opts.s_sot_min, opts.s_sot_max, opts.s_sot0};
            end

            % TODO(@anton) This _severely_ hurts performance over the vectorized assignment by doing N_stages vertcats of
            %              casadi symbolics vs just a vectorized assignment which does one. As such there needs to be backend
            %              work done for vdx to cache vertcats of SX somehow. Current theory is one can simply keep a queue of
            %              symbolics to be added in a cell array until a read is done, at which point we call a single vertcat
            %              on the whole queue which is _significantly_ faster.
            % 2d vars: Variables that are defined for each finite element.
            for ii=1:opts.N_stages
                % other derived values
                h0 = opts.h_k(ii);
                if obj.opts.use_fesd
                    ubh = (1 + opts.gamma_h) * h0;
                    lbh = (1 - opts.gamma_h) * h0;
                    if opts.time_rescaling && ~opts.use_speed_of_time_variables
                        % if only time_rescaling is true, speed of time and step size all lumped together, e.g., \hat{h}_{k,i} = s_n * h_{k,i}, hence the bounds need to be extended.
                        ubh = (1+opts.gamma_h)*h0*opts.s_sot_max;
                        lbh = (1-opts.gamma_h)*h0/opts.s_sot_min;
                    end
                    obj.w.h(ii,1:opts.N_finite_elements(ii)) = {{'h', 1}, lbh, ubh, h0};
                end
                if obj.opts.step_equilibration == StepEquilibrationMode.linear_complementarity
                    obj.w.B_max(ii,2:opts.N_finite_elements(ii)) = {{'B_max', dims.n_lambda},-inf,inf};
                    % forward sum of c(x) or backward sum of c(x).
                    obj.w.pi_c(ii,2:opts.N_finite_elements(ii)) = {{'pi_c', dims.n_c},-inf,inf};
                    % forward sum of lambda or backward sum of lambda.
                    obj.w.pi_lambda(ii,2:opts.N_finite_elements(ii)) = {{'pi_lambda', dims.n_lambda},-inf,inf};
                    % multipliers for c(x) in max kkt conditions.
                    obj.w.c_mult(ii,2:opts.N_finite_elements(ii)) = {{'c_mult', dims.n_c},0,inf};
                    % multipliers for lambda in max kkt conditions.
                    obj.w.lambda_mult(ii,2:opts.N_finite_elements(ii)) = {{'lambda_mult', dims.n_lambda},0,inf};
                    % pi_c and pi_lambda, i.e. vector of switch indicators.
                    obj.w.eta(ii,2:opts.N_finite_elements(ii)) = {{'eta', dims.n_lambda},0,inf};
                    % an or of all the switch indicators.
                    obj.w.nu(ii,2:opts.N_finite_elements(ii)) = {{'nu', 1},0,inf};
                end
            end

            % For c_n ~= 1 case
            rbp = ~opts.right_boundary_point_explicit;
            
            % 3d vars: Variables defined on each rk stage
            %          some of which are also defined at the initial point:
            obj.w.x(0,0,opts.n_s) = {{['x_0'], dims.n_x}, model.x0, model.x0, model.x0}; % differential state
            obj.w.z(0,0,opts.n_s) = {{'z', dims.n_z}, model.lbz, model.ubz, model.z0}; % user algebraics
            obj.w.lambda(0,0,opts.n_s) = {{['lambda'], dims.n_lambda},0,inf, opts.initial_lambda_gcs};
            obj.w.c_lift(0,0,opts.n_s) = {{['c_lift'], dims.n_c_lift},-inf,inf};
            for ii=1:opts.N_stages
                obj.w.x(ii,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp)) = {{'x', dims.n_x}, model.lbx, model.ubx, model.x0};
                if opts.rk_representation == RKRepresentation.differential_lift_x
                    obj.w.v(ii,1:opts.N_finite_elements(ii),1:opts.n_s) = {{'v', dims.n_x}}; % v = f_x
                end
                
                obj.w.z(ii,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp)) = {{'z', dims.n_z}, model.lbz, model.ubz, model.z0};
                obj.w.lambda(ii,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp)) = {{'lambda', dims.n_lambda},0, inf, opts.initial_lambda_gcs};
                obj.w.c_lift(ii,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp)) = {{'c_lift', dims.n_c_lift},-inf,inf};
                % Handle x_box settings
                if ~opts.x_box_at_stg
                    obj.w.x(1:opts.N_stages,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp-1)).lb = -inf*ones(dims.n_x, 1);
                    obj.w.x(1:opts.N_stages,1:opts.N_finite_elements(ii),1:(opts.n_s+rbp-1)).ub = inf*ones(dims.n_x, 1);
                end

                if ~opts.x_box_at_fe
                    obj.w.x(1:opts.N_stages,1:(opts.N_finite_elements(ii)-1),opts.n_s+rbp).lb = -inf*ones(dims.n_x, 1);
                    obj.w.x(1:opts.N_stages,1:(opts.N_finite_elements(ii)-1),opts.n_s+rbp).ub = inf*ones(dims.n_x, 1);
                end
            end
        end

        function generate_direct_transcription_constraints(obj)
            import casadi.*
            model = obj.model;
            opts = obj.opts;
            dcs = obj.dcs;
            dims = obj.dcs.dims;

            rbp = ~opts.right_boundary_point_explicit;

            v_global = obj.w.v_global();
            p_global = obj.p.p_global();

            % Recall that we treat (0,0,n_s) as the initial point. This is done for ergonomics in extracting results.
            x_0 = obj.w.x(0,0,opts.n_s);
            z_0 = obj.w.z(0,0,opts.n_s);
            c_lift_0 = obj.w.c_lift(0,0,opts.n_s);
            lambda_0 = obj.w.lambda(0,0,opts.n_s);
            
            obj.g.z(0,0,opts.n_s) = {dcs.g_z_fun(x_0, z_0, obj.w.u(1), v_global, [p_global;obj.p.p_time_var(1)])};
            obj.g.c_lb(0,0,opts.n_s) = {dcs.c_fun(x_0, c_lift_0, z_0, v_global, [p_global;obj.p.p_time_var(1)]), 0, inf};
            obj.g.c_lift(0,0,opts.n_s) = {dcs.g_c_lift_fun(x_0, c_lift_0, z_0, v_global, [p_global;obj.p.p_time_var(1)])};
            
            x_prev = obj.w.x(0,0,opts.n_s);
            for ii=1:opts.N_stages
                h0 = obj.p.T().val/(opts.N_stages*opts.N_finite_elements(ii));
                
                u_i = obj.w.u(ii);
                p_stage = obj.p.p_time_var(ii);
                p = [p_global;p_stage];
                if obj.opts.use_speed_of_time_variables && opts.local_speed_of_time_variable
                    s_sot = obj.w.sot(ii);
                elseif obj.opts.use_speed_of_time_variables
                    s_sot = obj.w.sot();
                else
                    s_sot = 1;
                end
                if opts.time_optimal_problem && ~opts.use_speed_of_time_variables
                    t_stage = obj.w.T_final()/(opts.N_stages*opts.N_finite_elements(ii));
                elseif opts.time_optimal_problem
                    t_stage = s_sot*obj.p.T()/opts.N_stages;
                else
                    t_stage = obj.p.T()/opts.N_stages;
                end

                sum_h = 0;
                for jj=1:opts.N_finite_elements(ii)
                    if obj.opts.use_fesd
                        h = obj.w.h(ii,jj);
                        sum_h = sum_h + h;
                    elseif opts.time_optimal_problem && ~opts.use_speed_of_time_variables
                        h = obj.w.T_final()/(opts.N_stages*opts.N_finite_elements(ii));
                    else
                        h = h0;
                    end
                    switch opts.rk_representation
                      case RKRepresentation.integral
                        % In integral representation stage variables are states.
                        x_ij_end = opts.D_rk(1)*x_prev;
                        for kk=1:opts.n_s
                            x_ijk = obj.w.x(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);
                            lambda_ijk = obj.w.lambda(ii,jj,kk);
                            c_lift_ijk = obj.w.c_lift(ii,jj,kk);

                            f_ijk = s_sot*dcs.f_x_fun(x_ijk, z_ijk, lambda_ijk, u_i, v_global, p);
                            q_ijk = s_sot*dcs.f_q_fun(x_ijk, z_ijk, lambda_ijk, u_i, v_global, p);
                            xk = opts.C_rk(1, kk+1) * x_prev;
                            for rr=1:opts.n_s
                                x_ijr = obj.w.x(ii,jj,rr);
                                xk = xk + opts.C_rk(rr+1, kk+1) * x_ijr;
                            end
                            obj.g.dynamics(ii,jj,kk) = {h * f_ijk - xk};
                            obj.g.z(ii,jj,kk) = {dcs.g_z_fun(x_ijk, z_ijk, u_i, v_global, p)};
                            obj.g.c_lb(ii,jj,kk) = {dcs.c_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p), 0, inf};
                            obj.g.c_lift(ii,jj,kk) = {dcs.g_c_lift_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p)};

                            x_ij_end = x_ij_end + opts.D_rk(kk+1)*x_ijk;
                            
                            if opts.g_path_at_stg
                                obj.g.path(ii,jj,kk) = {dcs.g_path_fun(x_ijk, z_ijk, u_i, v_global, p), model.lbg_path, model.ubg_path};
                            end
                            if size(model.G_path, 1) > 0
                                G_path = dcs.G_path_fun(x_ijk, z_ijk, u_i, v_global, p);
                                H_path = dcs.H_path_fun(x_ijk, z_ijk, u_i, v_global, p);
                                obj.G.path(ii,jj,kk) = {G_path};
                                obj.H.path(ii,jj,kk) = {H_path};
                            end
                            if ~opts.euler_cost_integration
                                % also integrate the objective
                                obj.f = obj.f + opts.B_rk(kk+1)*h*q_ijk;
                            end
                        end
                        if ~opts.right_boundary_point_explicit
                            error("not implemented")
                        end
                        if ~opts.g_path_at_stg && opts.g_path_at_fe
                            obj.g.path(ii,jj) = {dcs.g_path_fun(x_ijk, z_ijk, u_i, v_global, p), model.lbg_path, model.ubg_path};
                        end
                      case RKRepresentation.differential
                        error("Differential representation without lifting is unsupported for gradient complementarity systems")
                        
                      case RKRepresentation.differential_lift_x
                        % In differential representation with lifted state stage variables are the state derivatives and we
                        % lift the states at each stage point as well.
                        for kk = 1:opts.n_s
                            x_ijk = obj.w.x(ii,jj,kk);
                            x_temp = x_prev;
                            for rr = 1:opts.n_s
                                x_temp = x_temp + h*opts.A_rk(kk,rr)*obj.w.v(ii,jj,rr);
                            end
                            obj.g.lift_x(ii,jj,kk) = {x_ijk - x_temp};
                        end
                        x_ij_end = x_prev;
                        for kk=1:opts.n_s
                            x_ijk = obj.w.x(ii,jj,kk);
                            v_ijk = obj.w.v(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);
                            lambda_ijk = obj.w.lambda(ii,jj,kk);
                            c_lift_ijk = obj.w.c_lift(ii,jj,kk);
                           
                            f_ijk = s_sot*dcs.f_x_fun(x_ijk, z_ijk, lambda_ijk, u_i, v_global, p);
                            q_ijk = s_sot*dcs.f_q_fun(x_ijk, z_ijk, lambda_ijk, u_i, v_global, p);

                            x_ij_end = x_ij_end + h*opts.b_rk(kk)*v_ijk;
                            obj.g.v(ii,jj,kk) = {f_ijk - v_ijk};
                            obj.g.z(ii,jj,kk) = {dcs.g_z_fun(x_ijk, z_ijk, u_i, v_global, p)};
                            obj.g.c_lb(ii,jj,kk) = {dcs.c_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p), 0, inf};
                            obj.g.c_lift(ii,jj,kk) = {dcs.g_c_lift_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p)};
                            
                            if opts.g_path_at_stg
                                obj.g.path(ii,jj,kk) = {dcs.g_path_fun(x_ijk, z_ijk, u_i, v_global, p), model.lbg_path, model.ubg_path};
                            end
                            if size(model.G_path, 1) > 0
                                G_path = dcs.G_path_fun(x_ijk, z_ijk, u_i, v_global, p);
                                H_path = dcs.H_path_fun(x_ijk, z_ijk, u_i, v_global, p);
                                obj.G.path(ii,jj,kk) = {G_path};
                                obj.H.path(ii,jj,kk) = {H_path};
                            end
                            if ~opts.euler_cost_integration
                                % also integrate the objective
                                obj.f = obj.f + opts.b_rk(kk)*h*q_ijk;
                            end
                        end
                        if ~opts.g_path_at_stg && opts.g_path_at_fe
                            obj.g.path(ii,jj) = {dcs.g_path_fun(x_ijk, z_ijk, u_i, v_global, p), model.lbg_path, model.ubg_path};
                        end
                    end
                    x_prev = obj.w.x(ii,jj,opts.n_s+rbp);
                end
                if ~opts.g_path_at_stg && ~opts.g_path_at_fe
                    x_i = obj.w.x(ii, opts.N_finite_elements(ii), opts.n_s);
                    z_i = obj.w.z(ii, opts.N_finite_elements(ii), opts.n_s);
                    obj.g.path(ii) = {dcs.g_path_fun(x_i, z_i, u_i, v_global, p), model.lbg_path, model.ubg_path};
                end

                % Least Squares Costs
                % TODO we should convert the refs to params
                if ~isempty(model.x_ref_val)
                    obj.f = obj.f + t_stage*dcs.f_lsq_x_fun(obj.w.x(ii,opts.N_finite_elements(ii),opts.n_s),...
                        model.x_ref_val(:,ii),...
                        p);
                end
                if ~isempty(model.u_ref_val)
                    obj.f = obj.f + t_stage*dcs.f_lsq_u_fun(obj.w.u(ii),...
                        model.u_ref_val(:,ii),...
                        p);
                end
                if opts.euler_cost_integration
                    obj.f = obj.f + dcs.f_q_fun(x_ijk, z_ijk, lambda_ijk, u_i, v_global, p);
                end

                % Clock Constraints
                if obj.opts.use_fesd && opts.equidistant_control_grid
                    if opts.time_optimal_problem
                        if opts.use_speed_of_time_variables
                            obj.g.equidistant_control_grid(ii) = {[sum_h - opts.h;sot*sum_h - obj.w.T_final()/opts.N_stages]};
                        else
                            obj.g.equidistant_control_grid(ii) = {sum_h - obj.w.T_final()/opts.N_stages};
                        end
                    else
                        if opts.relax_terminal_numerical_time
                            % TODO(@armin) why is this allowed to be negative?
                            obj.w.s_numerical_time(ii) = {{'s_numerical', 1}, -2*opts.h, 2*opts.h, opts.h/2};
                            g_eq_grid = [sum_h - t_stage - obj.w.s_numerical_time(ii);
                                -(sum_h - t_stage) - obj.w.s_numerical_time(ii)];
                            obj.g.equidistant_control_grid(ii) = {g_eq_grid, -inf, 0};
                            obj.f = obj.f + opts.rho_terminal_numerical_time*obj.w.s_numerical_time(ii);
                        else
                            obj.g.equidistant_control_grid(ii) = {t_stage-sum_h};
                        end
                    end
                end
            end

            x_end = obj.w.x(opts.N_stages,opts.N_finite_elements(opts.N_stages),opts.n_s+rbp);
            z_end = obj.w.z(opts.N_stages,opts.N_finite_elements(opts.N_stages),opts.n_s+rbp);
            % Terminal cost
            obj.f = obj.f + dcs.f_q_T_fun(x_end, z_end, v_global, p_global);

            % Terminal_lsq_cost
            if ~isempty(model.x_ref_end_val)
                obj.f = obj.f + h0*opts.N_finite_elements(ii)*dcs.f_lsq_T_fun(x_end,...
                    model.x_ref_end_val,...
                    p_global);
            end

            % Terminal constraint
            if opts.relax_terminal_constraint_homotopy
                error("Currently unsupported")
            end
            g_terminal = dcs.g_terminal_fun(x_end, z_end, v_global, p_global);
            switch opts.relax_terminal_constraint
              case ConstraintRelaxationMode.NONE % hard constraint
                if opts.relax_terminal_constraint_from_above
                    obj.g.terminal = {g_terminal, model.lbg_terminal, inf*ones(dims.n_g_terminal,1)};
                else
                    obj.g.terminal = {g_terminal, model.lbg_terminal, model.ubg_terminal};
                end
              case ConstraintRelaxationMode.ELL_1 % l_1
                obj.w.s_terminal_ell_1 = {{'s_terminal_ell_1', dims.n_g_terminal}, 0, inf, 10};

                g_terminal = [g_terminal-model.lbg_terminal-obj.w.s_terminal_ell_1();
                    -(g_terminal-model.ubg_terminal)-obj.w.s_terminal_ell_1()];
                obj.g.terminal = {g_terminal, -inf, 0};
                obj.f = obj.f + obj.p.rho_terminal_p()*sum(obj.w.s_terminal_ell_1());
              case ConstraintRelaxationMode.ELL_2 % l_2
                     % TODO(@anton): this is as it was implemented before. should handle lb != ub?
                obj.f = obj.f + obj.p.rho_terminal_p()*(g_terminal-model.lbg_terminal)'*(g_terminal-model.lbg_terminal);
              case ConstraintRelaxationMode.ELL_INF % l_inf
                obj.w.s_terminal_ell_inf = {{'s_terminal_ell_inf', 1}, 0, inf, 1e3};

                g_terminal = [g_terminal-model.lbg_terminal-obj.w.s_terminal_ell_inf();
                    -(g_terminal-model.ubg_terminal)-obj.w.s_terminal_ell_inf()];
                obj.g.terminal = {g_terminal, -inf, 0};
                obj.f = obj.f + obj.p.rho_terminal_p()*obj.w.s_terminal_ell_inf();
            end
        end

        function generate_complementarity_constraints(obj)
            import casadi.*
            opts = obj.opts;
            dcs = obj.dcs;
            model = obj.model;
            % Do Cross-Complementarity

            rbp = ~opts.right_boundary_point_explicit;

            v_global = obj.w.v_global();
            p_global = obj.p.p_global();

            x_prev = obj.w.x(0,0,opts.n_s);
            z_prev = obj.w.z(0,0,opts.n_s);
            lambda_prev = obj.w.lambda(0,0,opts.n_s);
            c_lift_prev = obj.w.c_lift(0,0,opts.n_s);
            c_prev = dcs.c_fun(x_prev, c_lift_prev, z_prev, v_global, [p_global;obj.p.p_time_var(1)]);

            if opts.use_fesd
                switch opts.cross_comp_mode
                  case CrossCompMode.STAGE_STAGE
                    for ii=1:opts.N_stages
                        p_stage = obj.p.p_time_var(ii);
                        p = [p_global;p_stage];
                        for jj=1:opts.N_finite_elements(ii);
                            Gij = {};
                            Hij = {};
                            for rr=1:(opts.n_s + rbp)
                                x_ijr = obj.w.x(ii,jj,rr);
                                z_ijr = obj.w.z(ii,jj,rr);
                                c_lift_ijr = obj.w.c_lift(ii,jj,rr);
                                c_ijr = dcs.c_fun(x_ijr, c_lift_ijr, z_ijr, v_global, p);

                                Gij = vertcat(Gij, {lambda_prev});
                                Hij = vertcat(Hij, {c_ijr});
                            end
                            for rr=1:(opts.n_s + rbp)
                                lambda_ijr = obj.w.lambda(ii,jj,rr);

                                Gij = vertcat(Gij, {lambda_ijr});
                                Hij = vertcat(Hij, {c_prev});
                            end
                            for kk=1:(opts.n_s + rbp)
                                lambda_ijk = obj.w.lambda(ii,jj,kk);
                                for rr=1:(opts.n_s + rbp)
                                    x_ijr = obj.w.x(ii,jj,rr);
                                    z_ijr = obj.w.z(ii,jj,rr);
                                    c_lift_ijr = obj.w.c_lift(ii,jj,rr);
                                    c_ijr = dcs.c_fun(x_ijr, c_lift_ijr, z_ijr, v_global, p);
                                    
                                    Gij = vertcat(Gij, {lambda_ijk});
                                    Hij = vertcat(Hij, {c_ijr});
                                end
                            end
                            obj.G.cross_comp(ii,jj) = {vertcat(Gij{:})};
                            obj.H.cross_comp(ii,jj) = {vertcat(Hij{:})};
                            lambda_prev = obj.w.lambda(ii,jj,opts.n_s + rbp);
                            x_prev = obj.w.x(ii,jj,opts.n_s + rbp);
                            z_prev = obj.w.z(ii,jj,opts.n_s + rbp);
                            c_lift_prev = obj.w.c_lift(ii,jj,opts.n_s + rbp);
                            c_prev = dcs.c_fun(x_prev, c_lift_prev, z_prev, v_global, p);
                        end
                    end
                  case CrossCompMode.FE_STAGE
                    for ii=1:opts.N_stages
                        p_stage = obj.p.p_time_var(ii);
                        p = [p_global;p_stage];
                        for jj=1:opts.N_finite_elements(ii);
                            sum_lambda = lambda_prev + sum2(obj.w.lambda(ii,jj,:));
                            Gij = {sum_lambda};
                            Hij = {c_prev};
                            for kk=1:(opts.n_s + rbp)
                                x_ijk = obj.w.x(ii,jj,kk);
                                z_ijk = obj.w.z(ii,jj,kk);
                                c_lift_ijk = obj.w.c_lift(ii,jj,kk);
                                c_ijk = dcs.c_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p);

                                Gij = vertcat(Gij, {sum_lambda});
                                Hij = vertcat(Hij, {c_ijk});
                            end
                            obj.G.cross_comp(ii,jj) = {vertcat(Gij{:})};
                            obj.H.cross_comp(ii,jj) = {vertcat(Hij{:})};
                            lambda_prev = obj.w.lambda(ii,jj,opts.n_s + rbp);
                            x_prev = obj.w.x(ii,jj,opts.n_s + rbp);
                            z_prev = obj.w.z(ii,jj,opts.n_s + rbp);
                            c_lift_prev = obj.w.c_lift(ii,jj,opts.n_s + rbp);
                            c_prev = dcs.c_fun(x_prev, c_lift_prev, z_prev, v_global, p);
                        end
                    end
                  case CrossCompMode.STAGE_FE
                    for ii=1:opts.N_stages
                        p_stage = obj.p.p_time_var(ii);
                        p = [p_global;p_stage];
                        for jj=1:opts.N_finite_elements(ii);
                            sum_c = dcs.c_fun(x_prev);
                            for kk=1:(opts.n_s + rbp)
                                x_ijk = obj.w.x(ii,jj,kk);
                                z_ijk = obj.w.z(ii,jj,kk);
                                c_lift_ijk = obj.w.c_lift(ii,jj,kk);
                                c_ijk = dcs.c_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p);
                                sum_c = sum_c + c_ijk;
                            end
                            Gij = {lambda_prev};
                            Hij = {sum_c};
                            for kk=1:(opts.n_s + rbp)
                                lambda_ijk = obj.w.lambda(ii,jj,kk);

                                Gij = vertcat(Gij, {lambda_ijk});
                                Hij = vertcat(Hij, {sum_c});
                            end
                            obj.G.cross_comp(ii,jj) = {vertcat(Gij{:})};
                            obj.H.cross_comp(ii,jj) = {vertcat(Hij{:})};
                            lambda_prev = obj.w.lambda(ii,jj,opts.n_s + rbp);
                            x_prev = obj.w.x(ii,jj,opts.n_s + rbp);
                            z_prev = obj.w.z(ii,jj,opts.n_s + rbp);
                            c_lift_prev = obj.w.c_lift(ii,jj,opts.n_s + rbp);
                            c_prev = dcs.c_fun(x_prev, c_lift_prev, z_prev, v_global, p);
                        end
                    end
                  case CrossCompMode.FE_FE
                    for ii=1:opts.N_stages
                        p_stage = obj.p.p_time_var(ii);
                        p = [p_global;p_stage];
                        for jj=1:opts.N_finite_elements(ii);
                            sum_lambda = lambda_prev + sum2(obj.w.lambda(ii,jj,:));
                            sum_c = c_prev;
                            for kk=1:(opts.n_s + rbp)
                                x_ijk = obj.w.x(ii,jj,kk);
                                z_ijk = obj.w.z(ii,jj,kk);
                                c_lift_ijk = obj.w.c_lift(ii,jj,kk);
                                c_ijk = dcs.c_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p);
                                sum_c = sum_c + c_ijk;
                            end
                            obj.G.cross_comp(ii,jj) = {sum_lambda};
                            obj.H.cross_comp(ii,jj) = {sum_c};
                            lambda_prev = obj.w.lambda(ii,jj,opts.n_s + rbp);
                            x_prev = obj.w.x(ii,jj,opts.n_s + rbp);
                            z_prev = obj.w.z(ii,jj,opts.n_s + rbp);
                            c_lift_prev = obj.w.c_lift(ii,jj,opts.n_s + rbp);
                            c_prev = dcs.c_fun(x_prev, c_lift_prev, z_prev, v_global, p);
                        end
                    end
                end
            else
                obj.G.standard_comp(0,0,opts.n_s) = {lambda_prev};
                obj.H.standard_comp(0,0,opts.n_s) = {c_prev};
                for ii=1:opts.N_stages
                    p_stage = obj.p.p_time_var(ii);
                    p = [p_global;p_stage];
                    for jj=1:opts.N_finite_elements(ii);
                        Gij = {};
                        Hij = {};
                        for kk=1:(opts.n_s + rbp)
                            lambda_ijk = obj.w.lambda(ii,jj,kk);
                            x_ijk = obj.w.x(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);
                            c_lift_ijk = obj.w.c_lift(ii,jj,kk);
                            c_ijk = dcs.c_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p);

                            obj.G.standard_comp(ii,jj, kk) = {lambda_ijk};
                            obj.H.standard_comp(ii,jj, kk) = {c_ijk};
                        end
                    end
                end
            end
        end

        function generate_step_equilibration_constraints(obj)
            import casadi.*
            model = obj.model;
            dcs = obj.dcs;
            opts = obj.opts;
            dims = obj.dcs.dims;
            v_global = obj.w.v_global();
            p_global = obj.p.p_global();

            if ~opts.use_fesd % do nothing
                return
            end
            rbp = ~opts.right_boundary_point_explicit;

            switch obj.opts.step_equilibration
              case StepEquilibrationMode.heuristic_mean
                for ii=1:opts.N_stages
                    for jj=1:opts.N_finite_elements(ii)
                        h0 = obj.p.T()/(opts.N_stages*opts.N_finite_elements(ii));
                        obj.f = obj.f + obj.p.rho_h_p()*(h0-obj.w.h(ii,jj))^2;
                    end
                end
              case StepEquilibrationMode.heuristic_diff
                for ii=1:opts.N_stages
                    for jj=2:opts.N_finite_elements(ii)
                        h0 = obj.p.T()/(opts.N_stages*opts.N_finite_elements(ii));
                        obj.f = obj.f + obj.p.rho_h_p()*(obj.w.h(ii,jj)-obj.w.h(ii,jj-1))^2;
                    end
                end
              case StepEquilibrationMode.l2_relaxed_scaled
                eta_vec = [];
                for ii=1:opts.N_stages
                    p_stage = obj.p.p_time_var(ii);
                    p =[p_global;p_stage];
                    for jj=2:opts.N_finite_elements(ii)
                        sigma_c_B = 0;
                        sigma_lambda_B = 0;
                        for kk=1:(opts.n_s + rbp)
                            x_ijk = obj.w.x(ii,jj-1,kk);
                            z_ijk = obj.w.z(ii,jj-1,kk);
                            c_lift_ijk = obj.w.c_lift(ii,jj-1,kk);
                            c_ijk = dcs.c_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p);
                            sigma_c_B = sigma_c_B + c_ijk;
                            sigma_lambda_B = sigma_lambda_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lambda_F = 0;
                        for kk=1:(opts.n_s + rbp)
                            x_ijk = obj.w.x(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);
                            c_lift_ijk = obj.w.c_lift(ii,jj,kk);
                            c_ijk = dcs.c_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p);

                            sigma_c_F = sigma_c_F + c_ijk;
                            sigma_lambda_F = sigma_lambda_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lambda_B .* sigma_lambda_F;
                        nu = pi_c + pi_lam;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        obj.f = obj.f + obj.p.rho_h_p() * tanh(eta/opts.step_equilibration_sigma) * delta_h.^2;
                    end
                                end
              case StepEquilibrationMode.l2_relaxed
                eta_vec = [];
                for ii=1:opts.N_stages
                    p_stage = obj.p.p_time_var(ii);
                    p =[p_global;p_stage];
                    for jj=2:opts.N_finite_elements(ii)
                        sigma_c_B = 0;
                        sigma_lambda_B = 0;
                        for kk=1:(opts.n_s + rbp)
                            x_ijk = obj.w.x(ii,jj-1,kk);
                            z_ijk = obj.w.z(ii,jj-1,kk);
                            c_lift_ijk = obj.w.c_lift(ii,jj-1,kk);
                            c_ijk = dcs.c_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p);
                            sigma_c_B = sigma_c_B + c_ijk;
                            sigma_lambda_B = sigma_lambda_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lambda_F = 0;
                        for kk=1:(opts.n_s + rbp)
                            x_ijk = obj.w.x(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);
                            c_lift_ijk = obj.w.c_lift(ii,jj,kk);
                            c_ijk = dcs.c_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p);

                            sigma_c_F = sigma_c_F + c_ijk;
                            sigma_lambda_F = sigma_lambda_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lambda_B .* sigma_lambda_F;
                        nu = pi_c + pi_lam;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        obj.f = obj.f + obj.p.rho_h_p() * eta * delta_h.^2
                    end
                end
              case StepEquilibrationMode.direct
                eta_vec = [];
                for ii=1:opts.N_stages
                    p_stage = obj.p.p_time_var(ii);
                    p =[p_global;p_stage];
                    for jj=2:opts.N_finite_elements(ii)
                        sigma_c_B = 0;
                        sigma_lambda_B = 0;
                        for kk=1:(opts.n_s + rbp)
                            x_ijk = obj.w.x(ii,jj-1,kk);
                            z_ijk = obj.w.z(ii,jj-1,kk);
                            c_lift_ijk = obj.w.c_lift(ii,jj-1,kk);
                            c_ijk = dcs.c_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p);
                            sigma_c_B = sigma_c_B + c_ijk;
                            sigma_lambda_B = sigma_lambda_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lambda_F = 0;
                        for kk=1:(opts.n_s + rbp)
                            x_ijk = obj.w.x(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);
                            c_lift_ijk = obj.w.c_lift(ii,jj,kk);
                            c_ijk = dcs.c_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p);

                            sigma_c_F = sigma_c_F + c_ijk;
                            sigma_lambda_F = sigma_lambda_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lambda_B .* sigma_lambda_F;
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
                %obj.eta_fun = Function('eta_fun', {obj.w.sym}, {eta_vec});
              case StepEquilibrationMode.direct_homotopy
                error("not currently implemented")
                eta_vec = [];
                for ii=1:opts.N_stages
                    p_stage = obj.p.p_time_var(ii);
                    p =[p_global;p_stage];
                    for jj=2:opts.N_finite_elements(ii)
                        sigma_c_B = 0;
                        sigma_lambda_B = 0;
                        for kk=1:(opts.n_s + rbp)
                            x_ijk = obj.w.x(ii,jj-1,kk);
                            z_ijk = obj.w.z(ii,jj-1,kk);
                            c_lift_ijk = obj.w.c_lift(ii,jj-1,kk);
                            c_ijk = dcs.c_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p);
                            sigma_c_B = sigma_c_B + c_ijk;
                            sigma_lambda_B = sigma_lambda_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lambda_F = 0;
                        for kk=1:(opts.n_s + rbp)
                            x_ijk = obj.w.x(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);
                            c_lift_ijk = obj.w.c_lift(ii,jj,kk);
                            c_ijk = dcs.c_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p);

                            sigma_c_F = sigma_c_F + c_ijk;
                            sigma_lambda_F = sigma_lambda_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lambda_B .* sigma_lambda_F;
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
                %obj.eta_fun = Function('eta_fun', {obj.w.sym}, {eta_vec});
              case StepEquilibrationMode.linear_complementarity
                for ii=1:opts.N_stages
                    p_stage = obj.p.p_time_var(ii);
                    p =[p_global;p_stage];
                    h0 = obj.p.T()/(opts.N_stages*opts.N_finite_elements(ii));
                    for jj=2:opts.N_finite_elements(ii)
                        sigma_c_B = 0;
                        sigma_lambda_B = 0;
                        for kk=1:(opts.n_s + rbp)
                            x_ijk = obj.w.x(ii,jj-1,kk);
                            z_ijk = obj.w.z(ii,jj-1,kk);
                            c_lift_ijk = obj.w.c_lift(ii,jj-1,kk);
                            c_ijk = dcs.c_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p);
                            sigma_c_B = sigma_c_B + c_ijk;
                            sigma_lambda_B = sigma_lambda_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lambda_F = 0;
                        for kk=1:(opts.n_s + rbp)
                            x_ijk = obj.w.x(ii,jj,kk);
                            z_ijk = obj.w.z(ii,jj,kk);
                            c_lift_ijk = obj.w.c_lift(ii,jj,kk);
                            c_ijk = dcs.c_fun(x_ijk, c_lift_ijk, z_ijk, v_global, p);

                            sigma_c_F = sigma_c_F + c_ijk;
                            sigma_lambda_F = sigma_lambda_F + obj.w.lambda(ii,jj,kk);
                        end

                        lambda_mult = obj.w.lambda_mult(ii,jj);
                        c_mult = obj.w.c_mult(ii,jj);
                        B_max = obj.w.B_max(ii,jj);
                        pi_lambda = obj.w.pi_lambda(ii,jj);
                        pi_c = obj.w.pi_c(ii,jj);
                        eta = obj.w.eta(ii,jj);
                        nu = obj.w.nu(ii,jj);

                        obj.g.pi_lambda_or(ii,jj) = {[pi_lambda-sigma_lambda_F;pi_lambda-sigma_lambda_B;sigma_lambda_F+sigma_lambda_B-pi_lambda],0,inf};
                        obj.g.pi_c_or(ii,jj) = {[pi_c-sigma_c_F;pi_c-sigma_c_B;sigma_c_F+sigma_c_B-pi_c],0,inf};

                        % kkt conditions for min B, B>=sigmaB, B>=sigmaF
                        kkt_max = [1-c_mult-lambda_mult;
                            B_max-pi_lambda;
                            B_max-pi_c];
                        obj.g.kkt_max(ii,jj) = {kkt_max,
                            [0*ones(dims.n_lambda,1);0*ones(dims.n_lambda,1);0*ones(dims.n_lambda,1)],
                            [0*ones(dims.n_lambda,1);inf*ones(dims.n_lambda,1);inf*ones(dims.n_lambda,1)]};

                        obj.G.step_eq_kkt_max(ii,jj) = {[(B_max-pi_lambda);(B_max-pi_c)]};
                        obj.H.step_eq_kkt_max(ii,jj) = {[lambda_mult;c_mult]};
                        
                        % eta calculation
                        eta_const = [eta-pi_c;eta-pi_lambda;eta-pi_c-pi_lambda+B_max];
                        obj.g.eta_const(ii,jj) = {eta_const,
                            [-inf*ones(dims.n_lambda,1);-inf*ones(dims.n_lambda,1);zeros(dims.n_lambda,1)],
                            [zeros(dims.n_lambda,1);zeros(dims.n_lambda,1);inf*ones(dims.n_lambda,1)]};

                        obj.g.nu_or(ii,jj) = {[nu-eta;sum(eta)-nu],0,inf};

                        % the actual step eq conditions
                        M=opts.linear_complemtarity_M;
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        step_equilibration = [delta_h + (1/h0)*nu*M;
                            delta_h - (1/h0)*nu*M];
                        obj.g.step_equilibration(ii,jj) = {step_equilibration,[0;-inf],[inf;0]};
                    end
                end
            end
        end

        function populate_problem(obj)
            obj.create_variables();
            obj.generate_direct_transcription_constraints();
            obj.generate_complementarity_constraints();
            obj.generate_step_equilibration_constraints();

            obj.populated = true;
        end

        function create_solver(obj, solver_options, plugin, regenerate)
            if ~obj.populated
                obj.populate_problem();
            end

            if ~exist('plugin')
                plugin = 'scholtes_ineq';
            end

            if ~exist('regenerate')
                regenerate = false;
            end

            obj.finalize_assignments();
            
            % Sort by indices to recover almost block-band structure.
            % TODO: Maybe it would be useful to do a custom sorting for FESD to make the problem maximally block band.
            %       This is almost certainly necessary if we want to take advantage of e.g. FATROP. 
            if ~obj.sorted
                obj.w.sort_by_index();
                obj.g.sort_by_index();
            end
            solver_options.assume_lower_bounds = true;

            if regenerate || isempty(obj.solver) || (nosnoc.solver.MpccMethod(plugin) ~= obj.solver.relaxation_type)
                obj.solver = nosnoc.solver.mpccsol('Mpcc solver', plugin, obj.to_casadi_struct(), solver_options);
            end
        end

        function stats = solve(obj)
            opts = obj.opts;
            T_val = obj.p.T().val;

            if opts.use_fesd
                for ii=1:opts.N_stages
                    for jj=1:opts.N_finite_elements(ii)
                        % Recalculate ubh and lbh based on T_val
                        h0 = T_val/(opts.N_stages*opts.N_finite_elements(ii));
                        ubh = (1 + opts.gamma_h) * h0;
                        lbh = (1 - opts.gamma_h) * h0;
                        if opts.time_rescaling && ~opts.use_speed_of_time_variables
                            % if only time_rescaling is true, speed of time and step size all lumped together, e.g., \hat{h}_{k,i} = s_n * h_{k,i}, hence the bounds need to be extended.
                            ubh = (1+opts.gamma_h)*h0*opts.s_sot_max;
                            lbh = (1-opts.gamma_h)*h0/opts.s_sot_min;
                        end
                        obj.w.h(ii,jj).lb = lbh;
                        obj.w.h(ii,jj).ub = ubh;
                    end
                end
            end

            stats = solve@vdx.problems.Mpcc(obj);
        end
    end
end
