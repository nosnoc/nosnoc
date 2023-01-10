classdef NosnocProblem < NosnocFormulationObject
    properties
        % Index vectors
        ind_x
        ind_x0
        ind_u
        ind_v
        ind_theta
        ind_lam
        ind_mu
        ind_alpha
        ind_lambda_n
        ind_lambda_p
        ind_h
        ind_elastic
        ind_sot % index for speed of time variable
        ind_boundary % index of bundary value lambda and mu
        ind_t_final % Time-optimal problems: define auxilairy variable for the final time.

        ind_p_global
        ind_p_time_var

        ind_s_terminal

        % Problem data
        model
        settings
        dims
        ocp

        % Algorithmic parameters
        sigma_p
        rho_h_p
        rho_sot_p

        % Algorithmic global variables (time independent)
        s_elastic
        T_final

        % Parameters
        p
        p0

        % Problem components
        fe0
        stages % control stages

        % complementarity residual functions
        comp_res
        comp_std
        comp_fesd

        % Problem cost function
        cost_fun
    end
    % remaining list of TODOs
    % TODO: cleanup/add properties (in all components)
    % TODO: Create solver object, which will interact with setting parameters. 

    properties(Dependent, SetAccess=private, Hidden)
        u
        sot
        nu_vector
        
        ind_z
    end
    
    methods
        function obj = NosnocProblem(settings, dims, model)
            import casadi.*
            obj@NosnocFormulationObject();

            obj.ind_u = [];
            obj.ind_x = {};
            obj.ind_x0 = [];
            obj.ind_v = {};
            obj.ind_theta = {};
            obj.ind_lam = {};
            obj.ind_mu = {};
            obj.ind_alpha = {};
            obj.ind_lambda_n = {};
            obj.ind_lambda_p = {};
            obj.ind_h = {};
            obj.ind_sot = {};

            obj.ind_s_terminal = [];

            obj.settings = settings;
            obj.model = model;
            obj.dims = dims;
            obj.ocp = []; % TODO create ocp objects

            obj.stages = [];

            sigma_p = define_casadi_symbolic(settings.casadi_symbolic_mode, 'sigma_p');
            obj.sigma_p = sigma_p;
            rho_sot_p = define_casadi_symbolic(settings.casadi_symbolic_mode, 'rho_sot_p');
            obj.rho_sot_p = rho_sot_p;
            rho_h_p = define_casadi_symbolic(settings.casadi_symbolic_mode, 'rho_h_p');
            obj.rho_h_p = rho_h_p;
            rho_terminal_p = define_casadi_symbolic(settings.casadi_symbolic_mode, 'rho_terminal_p');
            T_ctrl_p  = define_casadi_symbolic(settings.casadi_symbolic_mode, 'T_ctrl_p');
            obj.p = [sigma_p;rho_sot_p;rho_h_p;rho_terminal_p;T_ctrl_p];

            % Populate parameter indices
            if dims.n_p_global > 0
                n_p = length(obj.p);
                ind_p_global = n_p+1:n_p+dims.n_p_global;
                obj.p = [obj.p; model.p_global];
            end
            if dims.n_p_time_var > 0;
                n_p = length(obj.p);
                ind_p_time_var = arrayfun(@(s) (n_p+(s*dims.n_p_time_var)+1):(n_p+(s*dims.n_p_time_var)+dims.n_p_time_var) , 0:dims.N_stages-1);
                obj.p = [obj.p; model.p_time_var];
            end

            if settings.time_optimal_problem
                % the final time in time optimal control problems
                T_final = define_casadi_symbolic(settings.casadi_symbolic_mode, 'T_final', 1);
                obj.T_final = T_final;
                T_final_guess = model.T;
            end

            % TODO Rename
            obj.createPrimalVariables();

            obj.createComplementarityConstraints();
            
            last_stage = obj.stages(end);
            last_fe = last_stage.stage(end);
            
            %  -- Constraint for the terminal numerical and physical time (if no equidistant grids are required) --
            % If the control grid is not equidistant, the constraint on sum of h happen only at the end.
            % The constraints are splited to those which are related to numerical and physical time, to make it easier to read.

            % terminal numerical and physical time
            % TODO: clean this up (sum_h/intergal_clock_state need to be functions probabaly)
            if settings.time_freezing
                % Terminal Phyisical Time (Possible terminal constraint on the clock state if time freezing is active).
                if settings.time_optimal_problem
                    obj.addConstraint(last_fe.x{end}(end)-T_final);
                else
                    if settings.impose_terminal_phyisical_time && ~settings.stagewise_clock_constraint
                        obj.addConstraint(last_fe.x{end}(end)-T_ctrl_p);
                    else
                        % no terminal constraint on the numerical time
                    end
                end
                if settings.equidistant_control_grid && ~settings.stagewise_clock_constraint
                    if ~settings.time_optimal_problem
                        obj.addConstraint(last_fe.x{end}(end)-model.T);
                    end
                end
            else
                if ~settings.use_fesd
                    if settings.time_optimal_problem
                        % if time_freezing is on, everything is done via the clock state.
                        if settings.use_speed_of_time_variables
                            integral_clock_state = 0;
                            for k=1:dims.N_stages
                                stage = obj.stages(k);
                                if settings.local_speed_of_time_variable
                                    s_sot = obj.sot{k};
                                else
                                    s_sot = obj.sot{1};
                                end
                                for fe=stage.stage
                                    integral_clock_state = integral_clock_state + fe.h*s_sot;
                                end
                            end
                            obj.addConstraint(integral_clock_state-T_final, 0, 0);
                        else
                            % otherwise treated via variable h_ki, i.e.,  h_ki =  T_final/(N_stages*N_FE)
                        end
                    end
                else
                    % if equidistant_control_grid = true all time constraint are added in
                    % the main control loop for every control stage k and the code
                    % below is skipped
                    if  ~settings.equidistant_control_grid
                        % T_num = T_phy = T_final =  T.
                        % all step sizes add up to prescribed time T.
                        % if use_speed_of_time_variables = true, numerical time is decupled from the sot scaling (no mather if local or not):
                        sum_h_all = 0;
                        for k=1:dims.N_stages
                            stage=obj.stages(k);
                            for fe=stage.stage
                                sum_h_all = sum_h_all+fe.h;
                            end
                        end
                        if ~settings.time_optimal_problem
                            obj.addConstraint(sum_h_all-model.T, 0, 0);
                        else
                            if ~settings.use_speed_of_time_variables
                                obj.addConstraint(sum_h_all-T_final, 0, 0);
                            else
                                integral_clock_state = 0;
                                for k=1:dims.N_stages
                                    stage = obj.stages(k);
                                    if settings.local_speed_of_time_variable
                                        s_sot = obj.sot{k};
                                    else
                                        s_sot = obj.sot{1};
                                    end
                                    for fe=stage.stage
                                        integral_clock_state = integral_clock_state + fe.h*s_sot;
                                    end
                                end
                                % T_num = T_phy = T_final \neq T.
                                obj.addConstraint(sum_h_all-model.T, 0, 0);
                                obj.addConstraint(integral_clock_state-T_final, 0, 0);
                            end
                        end
                    end
                end
            end

            % Process terminal constraint.
            if settings.terminal_constraint
                if settings.relax_terminal_constraint_homotopy
                    rho_terminal_p = 1/sigma_p;
                end
                X_end = last_fe.x{end};
                g_terminal = model.g_terminal_fun(X_end);
                n_terminal = length(g_terminal);
                if ~isequal(model.g_terminal_lb,model.g_terminal_ub)
                    settings.relax_terminal_constraint = 0;
                    if settings.print_level >2
                        fprintf('Info: Only terminal-equality constraint relaxation is supported, you have an inequality constraint.\n')
                    end
                end
                switch settings.relax_terminal_constraint % TODO name these.
                  case 0 % hard constraint
                    if settings.relax_terminal_constraint_from_above
                        obj.addConstraint(g_terminal, model.g_terminal_lb, inf*ones(n_terminal,1));
                    else
                        obj.addConstraint(g_terminal, model.g_terminal_lb, model.g_terminal_ub);
                    end
                  case 1 % l_1
                    s_terminal_ell_1 = define_casadi_symbolic(settings.casadi_symbolic_mode, 's_terminal_ell_1', n_terminal);
                    obj.addVariable(s_terminal_ell_1,...
                                    's_terminal',...
                                    1e3*ones(n_terminal,1),...
                                    -inf*ones(n_terminal,1),...
                                    inf*ones(n_terminal,1));

                    obj.addConstraint(g_terminal-g_terminal_lb-s_terminal_ell_1,...
                                      -inf*ones(n_terminal,1),...
                                      zeros(n_terminal,1));
                    obj.addConstraint(-(g_terminal-g_terminal_lb)-s_terminal_ell_1,...
                                      -inf*ones(n_terminal,1),...
                                      zeros(n_terminal,1));

                    obj.cost = obj.cost + rho_terminal_p*sum(s_terminal_ell_1);
                  case 2 % l_2
                    obj.cost = obj.cost + rho_terminal_p*(g_terminal-model.g_terminal_lb)'*(g_terminal-g_terminal_lb);
                  case 3 % l_inf
                    s_terminal_ell_inf = define_casadi_symbolic(settings.casadi_symbolic_mode, 's_terminal_ell_inf', 1);
                    obj.addVariable(s_terminal_ell_inf,...
                                    's_terminal',...
                                    1e3,...
                                    -inf,...
                                    inf);

                    obj.addConstraint(g_terminal-g_terminal_lb-s_terminal_ell_inf*ones(n_terminal,1),...
                                      -inf*ones(n_terminal,1),...
                                      zeros(n_terminal,1));
                    obj.addConstraint(-(g_terminal-g_terminal_lb)-s_terminal_ell_inf*ones(n_terminal,1),...
                                      -inf*ones(n_terminal,1),...
                                      zeros(n_terminal,1));

                    obj.cost = obj.cost + rho_terminal_p*s_terminal_ell_inf;
                  case 4 % l_inf, relaxed
                    % TODO: ask armin if this is correct.
                    if ismember(settings.mpcc_mode, MpccMode.elastic)
                        elastic = s_elastic*ones(n_terminal,1);
                    elseif ismemeber(settings.mpcc_mode, MpccMode.elastic_ell_1)
                        elastic = last_fe.elastic{end};
                    else
                        error('This mode of terminal constraint relaxation is only available if a MPCC elastic mode is used.');
                    end
                    obj.addConstraint(g_terminal-g_terminal_lb-elastic,...
                                      -inf*ones(n_terminal,1),...
                                      zeros(n_terminal,1));
                    obj.addConstraint(-(g_terminal-g_terminal_lb)-elastic,...
                                      -inf*ones(n_terminal,1),...
                                      zeros(n_terminal,1));
                end
            end

            % terminal least squares
            obj.cost = obj.cost + model.f_lsq_T_fun(last_fe.x{end}, model.x_ref_end_val);
            
            % Process terminal costs
            try
                obj.cost = obj.cost + model.f_q_T_fun(last_fe.x{end});
            catch
                fprintf('Terminal cost not defined');
            end
            
            % Process elastic costs
            if ismember(settings.mpcc_mode, MpccMode.elastic)
                if settings.objective_scaling_direct
                    obj.cost = obj.cost + (1/sigma_p)*obj.s_elastic;
                else
                    obj.cost = sigma_p*obj.cost + obj.s_elastic;
                end
            end
            if ismember(settings.mpcc_mode, MpccMode.elastic_ell_1)
                sum_s_elastic = 0;
                for k=1:dims.N_stages
                    stage=obj.stages(ii);
                    for fe=stage.stage
                        sum_s_elastic = sum_s_elastic + fe.sumElastic;
                    end
                end
                if settings.objective_scaling_direct
                    obj.cost = obj.cost + (1/sigma_p)*sum_s_elastic;
                else
                    obj.cost = sigma_p*obj.cost + sum_s_elastic;
                end
            end

            if settings.time_optimal_problem
                % Add to the vector of unknowns
                obj.addVariable(T_final, 't_final', settings.T_final_min, settings.T_final_max, T_final_guess);
                obj.cost = obj.cost + T_final;
            end

            % Calculate standard complementarities.
            J_comp_std = 0;
            for k=1:dims.N_stages
                stage = obj.stages(k);
                for fe=stage.stage
                    for j=1:dims.n_s
                        J_comp_std = J_comp_std + model.J_cc_fun(fe.rkStageZ(j));
                    end
                end
            end
            
            
            % Scalar-valued complementairity residual
            if settings.use_fesd
                J_comp_fesd = 0;
                for k=1:dims.N_stages
                    for fe=stage.stage
                        J_comp_fesd = J_comp_fesd + sum(diag(fe.sumTheta())*fe.sumLambda());
                    end
                end
                J_comp = J_comp_fesd;
            else
                J_comp_fesd = J_comp_std;
                J_comp = J_comp_std;
            end

            obj.comp_res = Function('comp_res', {obj.w, obj.p}, {J_comp});
            obj.comp_std = Function('comp_std', {obj.w, obj.p}, {J_comp_std});
            obj.comp_fesd = Function('comp_fesd', {obj.w, obj.p}, {J_comp_fesd});
            obj.cost_fun = Function('cost_fun', {obj.w}, {obj.cost});

            obj.p0 = [settings.sigma_0; settings.rho_sot; settings.rho_h; settings.rho_terminal; model.T];
            obj.p0

            if dims.n_p_global > 0;
                obj.p0 = [obj.p0; model.p_global_val];
            end
            
            if dims.n_p_time_var > 0;
                obj.p0 = [obj.p0; model.time_var_val];
            end
        end

        % TODO this should be private
        function createPrimalVariables(obj)
            import casadi.*
            fe0 = FiniteElementZero(obj.settings, obj.dims, obj.model);
            obj.fe0 = fe0;
            
            obj.p = vertcat(obj.p, fe0.lambda{1,:});

            obj.addVariable(fe0.x{1},...
                            'x0',...
                            fe0.lbw(fe0.ind_x{1}),...
                            fe0.ubw(fe0.ind_x{1}),...
                            fe0.w0(fe0.ind_x{1}));

            prev_fe = fe0;

            s_sot = [];
            if obj.settings.time_rescaling && obj.settings.use_speed_of_time_variables
                if ~obj.settings.local_speed_of_time_variable
                    s_sot = define_casadi_symbolic(obj.settings.casadi_symbolic_mode, 's_sot', 1);
                    obj.addVariable(s_sot,...
                                    'sot',...
                                    obj.settings.s_sot_min,...
                                    obj.settings.s_sot_max,...
                                    obj.settings.s_sot0);
                    if obj.settings.time_freezing
                        obj.cost = obj.cost + obj.rho_sot_p*(s_sot-1)^2;
                    end
                end
            end

            if ismember(obj.settings.mpcc_mode, MpccMode.elastic)
                s_elastic = define_casadi_symbolic(obj.settings.casadi_symbolic_mode, 's_elastic',1);
                obj.s_elastic = s_elastic;
                obj.addVariable(s_elastic, 'elastic', obj.settings.s_elastic_min, obj.settings.s_elastic_max, obj.settings.s_elastic_0);
            else
                s_elastic = [];
            end
                
            for ii=1:obj.dims.N_stages
                stage = ControlStage(prev_fe, obj.settings, obj.model, obj.dims, ii, s_sot, obj.T_final, obj.sigma_p, obj.rho_h_p, obj.rho_sot_p, s_elastic);
                obj.stages = [obj.stages, stage];

                obj.addControlStage(stage);
                obj.cost = obj.cost + stage.cost;
                prev_fe = stage.stage(end);
            end
        end

        function createComplementarityConstraints(obj)
            import casadi.*           
            model = obj.model;
            settings = obj.settings;
            dims = obj.dims;

            sigma_p = obj.sigma_p;
            s_elastic = obj.s_elastic;
            
            g_cross_comp = SX([]);
            % TODO Implement other modes
            if ~settings.use_fesd || settings.cross_comp_mode < 11
                % Do nothing, handled at the FE or stage level
                return
            elseif settings.cross_comp_mode == 11
                for r=1:dims.n_sys
                    g_r = 0;
                    for stage=obj.stages
                        for fe=stage.stage
                            g_r = g_r + diag(fe.sumTheta(r))*fe.sumLambda(r);
                        end
                    end
                    g_cross_comp = vertcat(g_cross_comp, g_r);
                end
            elseif settings.cross_comp_mode == 12
                for r=1:dims.n_sys
                    g_r = 0;
                    for stage=obj.stages
                        for fe=stage.stage
                            g_r = g_r + dot(fe.sumTheta(r),fe.sumLambda(r));
                        end
                    end
                    g_cross_comp = vertcat(g_cross_comp, g_r);
                end
            end

            g_comp = g_cross_comp;
            n_comp = length(g_cross_comp);
            
            %
            if ismember(settings.mpcc_mode, MpccMode.elastic_ell_1)
                s_elastic = define_casadi_symbolic(obj.settings.casadi_symbolic_mode, ['s_elastic_' num2str(obj.ctrl_idx) '_' num2str(obj.fe_idx)], n_comp);
                obj.addVariable(s_elastic,...
                                'elastic',...
                                settings.s_elastic_min*ones(n_comp,1),...
                                settings.s_elastic_max*ones(n_comp,1),...
                                settings.s_elastic_0*ones(n_comp,1));
            end
            
            % Do MPCC formulation
            % TODO This should be centralized
            if settings.mpcc_mode == MpccMode.Scholtes_ineq
                g_comp = g_comp - sigma_p;
                g_comp_ub = zeros(n_comp,1);
                g_comp_lb = -inf * ones(n_comp,1);
            elseif settings.mpcc_mode == MpccMode.Scholtes_eq
                g_comp = g_comp - sigma_p;
                g_comp_ub = zeros(n_comp,1);
                g_comp_lb = zeros(n_comp,1);
            elseif settings.mpcc_mode == MpccMode.ell_1_penalty
                if settings.objective_scaling_direct
                    obj.cost = obj.cost + (1/sigma_p)*sum(g_comp);
                else
                    obj.cost = sigma_p*obj.cost + sum(g_comp);
                end
            elseif settings.mpcc_mode == MpccMode.elastic_ineq
                g_comp = g_comp - s_elastic*ones(n_comp,1);
                g_comp_ub = zeros(n_comp,1);
                g_comp_lb = -inf * ones(n_comp,1);
            elseif settings.mpcc_mode == MpccMode.elastic_eq
                g_comp = g_comp - s_elastic*ones(n_comp,1);
                g_comp_ub = zeros(n_comp,1);
                g_comp_lb = zeros(n_comp,1);
            elseif settings.mpcc_mode == MpccMode.elastic_two_sided
                g_comp = [g_comp-s_elastic*ones(n_comp,1);g_comp+s_elastic*ones(n_comp,1)];
                g_comp_ub = [zeros(n_comp,1); inf*ones(n_comp,1)];
                g_comp_lb = [-inf*ones(n_comp,1);  zeros(n_comp,1)];
            elseif settings.mpcc_mode == MpccMode.elastic_ell_1_ineq
                g_comp = g_comp - s_elastic;
                g_comp_ub = zeros(n_comp,1);
                g_comp_lb = -inf * ones(n_comp,1);
            elseif settings.mpcc_mode == MpccMode.elastic_ell_1_eq
                g_comp = g_comp - s_elastic;
                g_comp_ub = zeros(n_comp,1);
                g_comp_lb = zeros(n_comp,1);
            elseif settings.mpcc_mode == MpccMode.elastic_ell_1_two_sided
                g_comp = [g_comp-s_elastic;g_comp+s_elastic];
                g_comp_ub = [zeros(n_comp,1); inf*ones(n_comp,1)];
                g_comp_lb = [-inf*ones(n_comp,1);  zeros(n_comp,1)];
            end

            if settings.mpcc_mode ~= MpccMode.ell_1_penalty
                obj.addConstraint(g_comp, g_comp_lb, g_comp_ub);
            end
        end
        
        % TODO this should be private
        function addControlStage(obj, stage)
            w_len = length(obj.w);

            obj.addPrimalVector(stage.w, stage.lbw, stage.ubw, stage.w0);

            obj.ind_h = [obj.ind_h, stage.ind_h+w_len];
            obj.ind_u = [obj.ind_u, stage.ind_u+w_len];
            obj.ind_sot = [obj.ind_sot, stage.ind_sot+w_len];
            obj.ind_x = [obj.ind_x, increment_indices(stage.ind_x, w_len)];
            obj.ind_v = [obj.ind_v, increment_indices(stage.ind_v, w_len)];
            obj.ind_theta = [obj.ind_theta, increment_indices(stage.ind_theta, w_len)];
            obj.ind_lam = [obj.ind_lam, increment_indices(stage.ind_lam, w_len)];
            obj.ind_mu = [obj.ind_mu, increment_indices(stage.ind_mu, w_len)];
            obj.ind_alpha = [obj.ind_alpha, increment_indices(stage.ind_alpha, w_len)];
            obj.ind_lambda_n = [obj.ind_lambda_n, increment_indices(stage.ind_lambda_n, w_len)];
            obj.ind_lambda_p = [obj.ind_lambda_p, increment_indices(stage.ind_lambda_p, w_len)];
            %obj.ind_nu_lift = [obj.ind_x, increment_indices(stage.ind_nu_lift, w_len)];

            obj.addConstraint(stage.g, stage.lbg, stage.ubg);
        end

        function addPrimalVector(obj, symbolic, lb, ub, initial)
            lens = [size(symbolic,1), size(lb,1), size(ub,1), size(initial,1)];
            if ~all(lens == lens(1))
                symbolic
                lb
                ub
                initial
                error("mismatched dims")
            end
            obj.w = vertcat(obj.w, symbolic);
            obj.lbw = vertcat(obj.lbw, lb);
            obj.ubw = vertcat(obj.ubw, ub);
            obj.w0 = vertcat(obj.w0, initial);
        end
        
        function u = get.u(obj)
            u = cellfun(@(u) obj.w(u), obj.ind_u, 'UniformOutput', false);
        end
        
        function sot = get.sot(obj)
            sot = cellfun(@(sot) obj.w(sot), obj.ind_sot, 'UniformOutput', false);
        end

        function nu_vector = get.nu_vector(obj)
            nu_vector = [];
            for k=obj.dims.N_stages
                stage = obj.stages(k);
                for fe=stage.stage
                    nu_vector = vertcat(nu_vector,fe.nu_vector);
                end
            end
        end

        function ind_z = get.ind_z(obj)
            % TODO Add gamma/beta
            ind_z = [flatten_ind(obj.ind_theta(1:obj.dims.n_s,:))
                     flatten_ind(obj.ind_lam(1:obj.dims.n_s,:))
                     flatten_ind(obj.ind_mu(1:obj.dims.n_s,:))
                     flatten_ind(obj.ind_alpha(1:obj.dims.n_s,:))
                     flatten_ind(obj.ind_lambda_n(1:obj.dims.n_s,:))
                     flatten_ind(obj.ind_lambda_p(1:obj.dims.n_s,:))];
            ind_z = sort(ind_z);
        end
        
        function print(obj)
            disp("g");
            disp(length(obj.g))
            print_casadi_vector(obj.g);
            disp('lbg, ubg');
            disp([length(obj.lbg), length(obj.ubg)]);
            disp([obj.lbg, obj.ubg]);

            disp("w");
            disp(length(obj.w))
            print_casadi_vector(obj.w);
            disp('lbw, ubw');
            disp([length(obj.lbw), length(obj.ubw)]);
            disp([obj.lbw, obj.ubw]);
            disp('w0');
            disp(length(obj.w0));
            disp(obj.w0);

            disp('objective');
            disp(obj.cost);
        end
    end
end
