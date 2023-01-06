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
        ind_nu_lift
        ind_h
        ind_elastic
        ind_sot % index for speed of time variable
        ind_boundary % index of bundary value lambda and mu
        ind_t_final % Time-optimal problems: define auxilairy variable for the final time.

        ind_s_terminal

        model
        settings
        dims
        ocp

        sigma_p

        rho_sot_p

        p
        p0

        fe0
        stages

        comp_res
        comp_std
        comp_fesd
        cost_fun
    end
    % remaining list of TODOs
    % TODO: cleanup properties

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

            obj.ind_u = cell(dims.N_stages,1);
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
            obj.ind_nu_lift = {};
            obj.ind_sot = {};

            obj.ind_s_terminal = [];

            obj.settings = settings;
            obj.model = model;
            obj.dims = dims;
            obj.ocp = []; % TODO create ocp objects

            obj.stages = {};

            sigma_p = SX.sym('sigma_p');
            rho_sot_p = SX.sym('rho_sot_p');
            obj.rho_sot_p = rho_sot_p;
            rho_h_p = SX.sym('rho_h_p');
            rho_terminal_p = SX.sym('rho_terminal_p');
            T_ctrl_p  = SX.sym('T_ctrl_p');
            obj.p = [sigma_p;rho_sot_p;rho_h_p;rho_terminal_p;T_ctrl_p];

            if settings.time_optimal_problem
                % the final time in time optimal control problems
                T_final = SX.sym('T_final', 1);
                T_final_guess = model.T;
            end

            if ismember(settings.mpcc_mode, MpccMode.elastic)
                s_elastic = SX.sym('s_elastic',1);
            else
                s_elastic = [];
            end
            obj.createPrimalVariables();

            for ii=1:dims.N_stages
                stage = obj.stages{ii};
                Uk = obj.u{ii};
                if settings.time_rescaling && settings.use_speed_of_time_variables
                    if settings.local_speed_of_time_variable
                        s_sot = obj.sot{ii};
                    else
                        s_sot = obj.sot{1};
                    end
                else
                    s_sot = 1;
                end

                obj.cost = obj.cost + (model.T/dims.N_stages)*model.f_lsq_x_fun(stage(end).x{end},model.x_ref_val(:,ii));
                if dims.n_u > 0
                    obj.cost = obj.cost + (model.T/dims.N_stages)*model.f_lsq_u_fun(Uk,model.u_ref_val(:,ii));
                end
                
                for fe=stage
                    % TODO: OCP
                    % 1) Stewart Runge-Kutta discretization
                    fe.forwardSimulation(obj.ocp, Uk, s_sot);

                    % 2) Complementarity Constraints
                    fe.createComplementarityConstraints(sigma_p, s_elastic);

                    % 3) Step Equilibration
                    fe.stepEquilibration(sigma_p, rho_h_p);

                    % 4) add finite element variables
                    obj.addFiniteElement(fe);
                    
                    % 5) add cost and constraints from FE to problem
                    obj.cost = obj.cost + fe.cost;
                    obj.addConstraint(fe.g, fe.lbg, fe.ubg);
                end

                % TODO: combine this into a function
                if settings.use_fesd && settings.equidistant_control_grid
                    if ~settings.time_optimal_problem
                        obj.addConstraint(sum(vertcat(stage.h)) - model.h);
                    elseif ~settings.time_freezing
                        if settings.use_speed_of_time_variables
                            obj.addConstraint(sum(vertcat(stage.h)) - model.h)
                            obj.addConstraint(sum(s_sot*vertcat(stage.h)) - T_final/dims.N_stages);
                        else
                            obj.addConstraint(sum(vertcat(stage.h)) - T_final/dims.N_stages);
                        end
                    end
                end
                if settings.time_freezing && settings.stagewise_clock_constraint
                    if time_optimal_problem
                        obj.addConstraint(fe.x{end}(end) - ii*(T_final/dims.N_stages) + model.x0(end));
                    else
                        obj.addConstraint(fe.x{end}(end) - ii*model.h + model.x0(end));
                    end
                end
            end
            last_fe = obj.stages{end}(end);
             % TODO time_freezing/time_optimal numerical/physical time

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
                    s_terminal_ell_1 = SX.sym('s_terminal_ell_1', n_terminal);
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
                    s_terminal_ell_inf = SX.sym('s_terminal_ell_inf', 1);
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

            %  -- Constraint for the terminal numerical and physical time (if no equidistant grids are required) --
            % If the control grid is not equidistant, the constraint on sum of h happen only at the end.
            % The constraints are splited to those which are related to numerical and physical time, to make it easier to read.

            % terminal numerical and physical time
            % TODO: clean this up (sum_h/intergal_clock_state need to be functions probabaly)
            if settings.time_freezing
                % Terminal Phyisical Time (Posssble terminal constraint on the clock state if time freezing is active).
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
                                stage = obj.stages{k};
                                if settings.local_speed_of_time_variable
                                    s_sot = obj.sot(ii);
                                else
                                    s_sot = obj.sot(1);
                                end
                                for fe=stage
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
                            stage=obj.stages{k};
                            for fe=stage
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
                                    stage = obj.stages{k};
                                    if settings.local_speed_of_time_variable
                                        s_sot = obj.sot(ii);
                                    else
                                        s_sot = obj.sot(1);
                                    end
                                    for fe=stage
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

            % terminal least squares
            obj.cost = obj.cost + model.f_lsq_T_fun(last_fe.x{end}, model.x_ref_end_val);
            
            % Process terminal costs
            try
                last_fe = obj.stages{end}(end);
                obj.cost = obj.cost + model.f_q_T_fun(last_fe.x(end));
            catch
                fprintf('Terminal cost not defined');
            end
            
            % Process elastic costs
            if ismember(settings.mpcc_mode, MpccMode.elastic)
                obj.addVariable(s_elastic, 'elastic', settings.s_elastic_min, settings.s_elastic_max, settings.s_elastic_0);
                if settings.objective_scaling_direct
                    obj.cost = obj.cost + (1/sigma_p)*s_elastic;
                else
                    obj.cost = sigma_p*obj.cost + s_elastic;
                end
            end
            if ismember(settings.mpcc_mode, MpccMode.elastic_ell_1)
                sum_s_elastic = 0;
                for k=1:dims.N_stages
                    stage=obj.stages{ii};
                    for fe=stage
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
                for fe=stage
                    for j=1:dims.n_s
                        J_comp_std = J_comp_std + model.J_cc_fun(fe.rkStageZ(j));
                    end
                end
            end
            
            
            % Scalar-valued complementairity residual
            if settings.use_fesd
                J_comp_fesd = 0;
                for k=1:dims.N_stages
                    for fe=stage
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

            obj.p0 = [settings.sigma_0, settings.rho_sot, settings.rho_h, settings.rho_terminal, model.T];
        end

        % TODO this should be private
        function createPrimalVariables(obj)
            fe0 = FiniteElementZero(obj.settings, obj.dims, obj.model);
            obj.fe0 = fe0;
            
            obj.p = vertcat(obj.p, fe0.lambda{1,:});

            obj.addVariable(fe0.x{1},...
                            'x0',...
                            fe0.lbw(fe0.ind_x{1}),...
                            fe0.ubw(fe0.ind_x{1}),...
                            fe0.w0(fe0.ind_x{1}));

            prev_fe = fe0;
            for ii=1:obj.dims.N_stages
                stage = obj.createControlStage(ii, prev_fe);
                obj.stages = [obj.stages, {stage}];
                prev_fe = stage(end);
            end
        end

        % TODO this should be private
        function control_stage = createControlStage(obj, ctrl_idx, prev_fe)
            import casadi.*
            Uk = SX.sym(['U_' num2str(ctrl_idx)], obj.dims.n_u);
            obj.addVariable(Uk, 'u', obj.model.lbu, obj.model.ubu,...
                            zeros(obj.dims.n_u,1), ctrl_idx);

            if obj.settings.time_rescaling && obj.settings.use_speed_of_time_variables
                if obj.settings.local_speed_of_time_variable
                    % at every stage
                    s_sot_k = SX.sym(['s_sot_' num2str(ctrl_idx)], 1);
                    obj.addVariable(s_sot_k,...
                                    'sot',...
                                    obj.settings.s_sot_min,...
                                    obj.settings.s_sot_max,...
                                    obj.settings.s_sot0,...
                                    ctrl_idx);
                    if obj.settings.time_freezing
                        obj.cost = obj.cost + obj.rho_sot_p*(s_sot_k-1)^2;
                    end
                else
                    if ctrl_idx == 1
                        % only once
                        s_sot = SX.sym('s_sot', 1);
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
            end
            
            control_stage = [];
             for ii=1:obj.dims.N_finite_elements
                fe = FiniteElement(prev_fe, obj.settings, obj.model, obj.dims, ctrl_idx, ii);
                control_stage = [control_stage, fe];
                prev_fe = fe;
            end
        end

        % TODO this should be private
        function addFiniteElement(obj, fe)
            w_len = length(obj.w);

            obj.addPrimalVector(fe.w, fe.lbw, fe.ubw, fe.w0);

            obj.ind_h = [obj.ind_h, fe.ind_h+w_len];
            obj.ind_x = [obj.ind_x, increment_indices(fe.ind_x, w_len)];
            obj.ind_v = [obj.ind_v, increment_indices(fe.ind_v, w_len)];
            obj.ind_theta = [obj.ind_theta, increment_indices(fe.ind_theta, w_len)];
            obj.ind_lam = [obj.ind_lam, increment_indices(fe.ind_lam, w_len)];
            obj.ind_mu = [obj.ind_mu, increment_indices(fe.ind_mu, w_len)];
            obj.ind_alpha = [obj.ind_alpha, increment_indices(fe.ind_alpha, w_len)];
            obj.ind_lambda_n = [obj.ind_lambda_n, increment_indices(fe.ind_lambda_n, w_len)];
            obj.ind_lambda_p = [obj.ind_lambda_p, increment_indices(fe.ind_lambda_p, w_len)];
            obj.ind_nu_lift = [obj.ind_x, fe.ind_nu_lift+w_len];
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
                stage = obj.stages{k};
                for fe=stage
                    nu_vector = vertcat(nu_vector,fe.nu_vector);
                end
            end
        end

        function ind_z = get.ind_z(obj)
            ind_z = [flatten_ind(obj.ind_theta)
                     flatten_ind(obj.ind_lam)
                     flatten_ind(obj.ind_mu)
                     flatten_ind(obj.ind_alpha)
                     flatten_ind(obj.ind_lambda_n)
                     flatten_ind(obj.ind_lambda_p)];
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

