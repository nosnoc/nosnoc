% BSD 2-Clause License

% Copyright (c) 2022, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% This file is part of NOSNOC.

classdef NosnocMPCC < NosnocFormulationObject
    properties % TODO: 
        % Index vectors
        ind_x
        ind_x0
        ind_u
        ind_v
        ind_z
        % Stewart
        ind_theta
        ind_lam
        ind_mu
        % Step
        ind_alpha
        ind_lambda_n
        ind_lambda_p
        ind_theta_step
        ind_beta
        % Speficif to CLS representation
        ind_lambda_normal
        ind_lambda_tangent
        ind_y_gap
        % friction multipliers and lifting
        % conic
        ind_gamma
        ind_beta_conic
        % poly
        ind_gamma_d
        ind_beta_d
        ind_delta_d
        % variables related to conic
        ind_p_vt
        ind_n_vt
        ind_alpha_vt
        % variables only at element boundary
        ind_x_left_bp
        ind_Y_gap
        ind_Lambda_normal
        ind_Lambda_tangent
        ind_Gamma
        ind_Gamma_d
        ind_Beta_conic
        ind_Beta_d
        ind_Delta_d
        ind_L_vn
        ind_P_vt
        ind_N_vt

        ind_Alpha_vt
        % misc
        ind_nu_lift
        ind_h
        ind_elastic
        ind_sot % index for speed of time variable
        ind_t_final % Time-optimal problems: define auxilairy variable for the final time.
        ind_v_global
        ind_s_terminal
        ind_comp_lift

        % Parameter index variables
        ind_p_x0
        ind_p_global
        ind_p_time_var

        % g at mpcc level
        ind_g_mpcc

        % Problem data
        model
        problem_options
        dims
        ocp

        % original initialization
        w0_original

        % Algorithmic parameters
        rho_h_p
        rho_sot_p

        % Algorithmic global variables (time independent)
        s_elastic
        T_final

        % Parameters
        p
        p0

        % cross comps
        cross_comps
        ind_std_comp

        % Problem components
        fe0 % Zeroth finite element (contains X0, lambda00)
        stages % control stages

        % complementarity residual functions
        comp_res
        comp_std
        comp_fesd

        % Problem cost function
        augmented_objective_fun

        % Problem objective function
        objective_fun

        % Problem constraint function
        g_fun

        %
        G_fun
        H_fun

        % switch indicator function
        nu_fun

        nabla_J
        nabla_J_fun
    end

    properties(Dependent, SetAccess=private, Hidden)
        % Properties generated on the fly.

        % casadi symbolics/expresions for u, sot, and nu
        u
        sot
        nu_vector
        cc_vector

        % Indices for all algebraic vars in the problem
        ind_z_all
        ind_x_all
    end

    methods
        function obj = NosnocMPCC(problem_options, model)
            import casadi.*
            tic;
            obj@NosnocFormulationObject();

            problem_options.preprocess();
            model.verify_and_backfill(problem_options);
            model.generate_variables(problem_options);
            model.generate_equations(problem_options);

            dims = model.dims;
            if problem_options.right_boundary_point_explicit || problem_options.dcs_mode == DcsMode.CLS
                rbp_allowance = 0;
            else
                rbp_allowance = 1;
            end

            if problem_options.dcs_mode == "CLS" && ~problem_options.right_boundary_point_explicit
                rbp_x_only = 1;
            else
                rbp_x_only = 0;
            end

            right_ygap = 0;
            if problem_options.dcs_mode == "CLS" && ~problem_options.right_boundary_point_explicit
                right_ygap = 1;
            end

            obj.ind_u = [];
            obj.ind_x0 = [];
            obj.ind_x = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance+rbp_x_only);
            obj.ind_v = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s);
            obj.ind_z = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            % Stewart
            obj.ind_theta = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_lam = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_mu = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            % Step
            obj.ind_alpha = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_lambda_n = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_lambda_p = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_beta = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_theta_step = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            % CLS
            obj.ind_lambda_normal = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_lambda_tangent = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_y_gap = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+right_ygap+rbp_allowance);
            % friction multipliers and lifting
            % conic
            obj.ind_gamma =  cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_beta_conic = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            % poly
            obj.ind_gamma_d = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_beta_d = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_delta_d = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            % variables related to conic
            obj.ind_p_vt = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_n_vt = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_alpha_vt = cell(problem_options.N_stages,problem_options.N_finite_elements(1),dims.n_s+rbp_allowance);

            % Impulse
            obj.ind_x_left_bp = cell(problem_options.N_stages,problem_options.N_finite_elements(1),1);
            obj.ind_Y_gap = cell(problem_options.N_stages,problem_options.N_finite_elements(1),1);
            obj.ind_Lambda_normal = cell(problem_options.N_stages,problem_options.N_finite_elements(1),1);
            obj.ind_Lambda_tangent = cell(problem_options.N_stages,problem_options.N_finite_elements(1),1);
            obj.ind_Gamma = cell(problem_options.N_stages,problem_options.N_finite_elements(1),1);
            obj.ind_Beta_conic = cell(problem_options.N_stages,problem_options.N_finite_elements(1),1);
            obj.ind_Gamma_d = cell(problem_options.N_stages,problem_options.N_finite_elements(1),1);
            obj.ind_Beta_d = cell(problem_options.N_stages,problem_options.N_finite_elements(1),1);
            obj.ind_Delta_d = cell(problem_options.N_stages,problem_options.N_finite_elements(1),1);
            obj.ind_L_vn = cell(problem_options.N_stages,problem_options.N_finite_elements(1),1);
            obj.ind_P_vt = cell(problem_options.N_stages,problem_options.N_finite_elements(1),1);
            obj.ind_N_vt = cell(problem_options.N_stages,problem_options.N_finite_elements(1),1);
            obj.ind_Alpha_vt = cell(problem_options.N_stages,problem_options.N_finite_elements(1),1);

            
            % misc
            obj.ind_nu_lift = {};
            obj.ind_h = {};
            obj.ind_sot = {};
            obj.ind_v_global = [];
            obj.ind_comp_lift = {};

            obj.ind_s_terminal = [];

            obj.problem_options = problem_options;
            obj.model = model;
            obj.dims = dims;
            obj.ocp = []; % TODO create ocp objects

            obj.stages = [];

            rho_sot_p = define_casadi_symbolic(problem_options.casadi_symbolic_mode, 'rho_sot_p');
            obj.rho_sot_p = rho_sot_p;
            rho_h_p = define_casadi_symbolic(problem_options.casadi_symbolic_mode, 'rho_h_p');
            obj.rho_h_p = rho_h_p;
            rho_terminal_p = define_casadi_symbolic(problem_options.casadi_symbolic_mode, 'rho_terminal_p');
            T_ctrl_p  = define_casadi_symbolic(problem_options.casadi_symbolic_mode, 'T_ctrl_p');
            obj.p = [rho_sot_p;rho_h_p;rho_terminal_p;T_ctrl_p];

            % Populate parameter indices
            if dims.n_p_global > 0
                n_p = length(obj.p);
                obj.ind_p_global = n_p+1:n_p+dims.n_p_global;
                obj.p = [obj.p; model.p_global];
            end
            if dims.n_p_time_var > 0;
                n_p = length(obj.p);
                obj.ind_p_time_var = arrayfun(@(s) (n_p+(s*dims.n_p_time_var)+1):(n_p+(s*dims.n_p_time_var)+dims.n_p_time_var) , 0:problem_options.N_stages-1);
                obj.p = [obj.p; model.p_time_var_stages(:)];
            end
            if problem_options.time_optimal_problem
                % the final time in time optimal control problems
                T_final = define_casadi_symbolic(problem_options.casadi_symbolic_mode, 'T_final', 1);
                obj.T_final = T_final;
                T_final_guess = problem_options.T;
            end

            % Add global vars
            obj.addVariable(model.v_global,...
                'v_global',...
                model.lbv_global,...
                model.ubv_global,...
                model.v0_global)

            obj.create_primal_variables();

            last_stage = obj.stages(end);
            last_fe = last_stage.stage(end);

            %  -- Constraint for the terminal numerical and physical time (if no equidistant grids are required) --
            % If the control grid is not equidistant, the constraint on sum of h happen only at the end.
            % The constraints are splited to those which are related to numerical and physical time, to make it easier to read.

            % terminal numerical and physical time
            % TODO: clean this up (sum_h/intergal_clock_state need to be functions probabaly)
            if problem_options.time_freezing
                % Terminal Phyisical Time (Possible terminal constraint on the clock state if time freezing is active).
                if problem_options.time_optimal_problem
                    obj.addConstraint(last_fe.x{end}(end)-T_final, 'type', 'g_mpcc');
                else
                    if problem_options.impose_terminal_phyisical_time && ~problem_options.stagewise_clock_constraint
                        obj.addConstraint(last_fe.x{end}(end)-T_ctrl_p, 'type', 'g_mpcc');
                    else
                        % no terminal constraint on the numerical time
                    end
                end
                if problem_options.equidistant_control_grid && ~problem_options.stagewise_clock_constraint
                    if ~problem_options.time_optimal_problem
                        obj.addConstraint(last_fe.x{end}(end)-problem_options.T, 'type', 'g_mpcc');
                    end
                end
            else
                if ~problem_options.use_fesd
                    if problem_options.time_optimal_problem
                        % if time_freezing is on, everything is done via the clock state.
                        if problem_options.use_speed_of_time_variables
                            integral_clock_state = 0;
                            for k=1:problem_options.N_stages
                                stage = obj.stages(k);
                                if problem_options.local_speed_of_time_variable
                                    s_sot = obj.sot{k};
                                else
                                    s_sot = obj.sot{1};
                                end
                                for fe=stage.stage
                                    integral_clock_state = integral_clock_state + fe.h*s_sot;
                                end
                            end
                            obj.addConstraint(integral_clock_state-T_final, 0, 0, 'type', 'g_mpcc');
                        else
                            % otherwise treated via variable h_ki, i.e.,  h_ki =  T_final/(N_stages*N_FE)
                        end
                    end
                else
                    % if equidistant_control_grid = true all time constraint are added in
                    % the main control loop for every control stage k and the code
                    % below is skipped
                    if  ~problem_options.equidistant_control_grid
                        % T_num = T_phy = T_final =  T.
                        % all step sizes add up to prescribed time T.
                        % if use_speed_of_time_variables = true, numerical time is decupled from the sot scaling (no mather if local or not):
                        sum_h_all = 0;
                        for k=1:problem_options.N_stages
                            stage=obj.stages(k);
                            for fe=stage.stage
                                sum_h_all = sum_h_all+fe.h;
                            end
                        end
                        if ~problem_options.time_optimal_problem
                            obj.addConstraint(sum_h_all-problem_options.T, 0, 0, 'type', 'g_mpcc');
                        else
                            if ~problem_options.use_speed_of_time_variables
                                obj.addConstraint(sum_h_all-T_final, 0, 0, 'type', 'g_mpcc');
                            else
                                integral_clock_state = 0;
                                for k=1:problem_options.N_stages
                                    stage = obj.stages(k);
                                    if problem_options.local_speed_of_time_variable
                                        s_sot = obj.sot{k};
                                    else
                                        s_sot = obj.sot{1};
                                    end
                                    for fe=stage.stage
                                        integral_clock_state = integral_clock_state + fe.h*s_sot;
                                    end
                                end
                                % T_num = T_phy = T_final \neq T.
                                obj.addConstraint(sum_h_all-problem_options.T, 0, 0, 'type', 'g_mpcc');
                                obj.addConstraint(integral_clock_state-T_final, 0, 0, 'type', 'g_mpcc');
                            end
                        end
                    end
                end
            end

            % Process terminal constraint.
            if ~isempty(model.g_terminal)
                if problem_options.relax_terminal_constraint_homotopy
                    rho_terminal_p = 1/sigma_p;
                end
                X_end = last_fe.x{end};
                g_terminal = model.g_terminal_fun(X_end, model.p_global, model.v_global);
                n_terminal = length(g_terminal);
                if ~isequal(model.g_terminal_lb,model.g_terminal_ub)
                    problem_options.relax_terminal_constraint = 0;
                    if problem_options.print_level >2
                        fprintf('Info: Only terminal-equality constraint relaxation is supported, you have an inequality constraint.\n')
                    end
                end
                switch problem_options.relax_terminal_constraint % TODO name these.
                  case 0 % hard constraint
                    if problem_options.relax_terminal_constraint_from_above
                        obj.addConstraint(g_terminal, model.g_terminal_lb, inf*ones(n_terminal,1), 'type', 'g_mpcc');
                    else
                        obj.addConstraint(g_terminal, model.g_terminal_lb, model.g_terminal_ub, 'type', 'g_mpcc');
                    end
                  case 1 % l_1
                    s_terminal_ell_1 = define_casadi_symbolic(problem_options.casadi_symbolic_mode, 's_terminal_ell_1', n_terminal);
                    obj.addVariable(s_terminal_ell_1,...
                        's_terminal',...
                        -inf*ones(n_terminal,1),...
                        inf*ones(n_terminal,1),...
                        1e3*ones(n_terminal,1));

                    obj.addConstraint(g_terminal-model.g_terminal_lb-s_terminal_ell_1,...
                        -inf*ones(n_terminal,1),...
                        zeros(n_terminal,1), 'type', 'g_mpcc');
                    obj.addConstraint(-(g_terminal-model.g_terminal_lb)-s_terminal_ell_1,...
                        -inf*ones(n_terminal,1),...
                        zeros(n_terminal,1), 'type', 'g_mpcc');

                    obj.augmented_objective = obj.augmented_objective + rho_terminal_p*sum(s_terminal_ell_1);
                  case 2 % l_2
                    obj.augmented_objective = obj.augmented_objective + rho_terminal_p*(g_terminal-model.g_terminal_lb)'*(g_terminal-model.g_terminal_lb);
                  case 3 % l_inf
                    s_terminal_ell_inf = define_casadi_symbolic(problem_options.casadi_symbolic_mode, 's_terminal_ell_inf', 1);
                    obj.addVariable(s_terminal_ell_inf,...
                        's_terminal',...
                        -inf,...
                        inf,...
                        1e3);

                    obj.addConstraint(g_terminal-model.g_terminal_lb-s_terminal_ell_inf*ones(n_terminal,1),...
                        -inf*ones(n_terminal,1),...
                        zeros(n_terminal,1), 'type', 'g_mpcc');
                    obj.addConstraint(-(g_terminal-model.g_terminal_lb)-s_terminal_ell_inf*ones(n_terminal,1),...
                        -inf*ones(n_terminal,1),...
                        zeros(n_terminal,1), 'type', 'g_mpcc');

                    obj.augmented_objective = obj.augmented_objective + rho_terminal_p*s_terminal_ell_inf;
                  case 4 % l_inf, relaxed
                    if ismember(problem_options.mpcc_mode, MpccMode.elastic)
                        elastic = s_elastic*ones(n_terminal,1);
                    elseif ismemeber(problem_options.mpcc_mode, MpccMode.elastic_ell_1)
                        elastic = last_fe.elastic{end};
                    else
                        error('This mode of terminal constraint relaxation is only available if a MPCC elastic mode is used.');
                    end
                    obj.addConstraint(g_terminal-model.g_terminal_lb-elastic,...
                        -inf*ones(n_terminal,1),...
                        zeros(n_terminal,1), 'type', 'g_mpcc');
                    obj.addConstraint(-(g_terminal-model.g_terminal_lb)-elastic,...
                        -inf*ones(n_terminal,1),...
                        zeros(n_terminal,1), 'type', 'g_mpcc');
                end
            end

            % terminal least squares
            obj.augmented_objective = obj.augmented_objective + model.f_lsq_T_fun(last_fe.x{end}, model.x_ref_end_val, model.p_global);
            obj.objective = obj.objective + model.f_lsq_T_fun(last_fe.x{end}, model.x_ref_end_val, model.p_global);

            % Process terminal costs
            obj.augmented_objective = obj.augmented_objective + model.f_q_T_fun(last_fe.x{end}, model.p_global, model.v_global);
            obj.objective = obj.objective + model.f_q_T_fun(last_fe.x{end}, model.p_global, model.v_global);

            if problem_options.time_optimal_problem
                % Add to the vector of unknowns
                obj.addVariable(T_final, 't_final', problem_options.T_final_min, problem_options.T_final_max, T_final_guess);
                obj.augmented_objective = obj.augmented_objective + T_final;
                obj.objective = obj.objective + T_final;
            end

            % Calculate standard complementarities.
            J_comp_std = 0;
            J_comp_std_infty = 0;
            for k=1:problem_options.N_stages
                stage = obj.stages(k);
                for fe=stage.stage
                    for j=1:dims.n_s
                        J_comp_std = J_comp_std + model.J_cc_fun(fe.rkStageZ(j));
                    end
                end
            end

            % calculate complementarity residual via vector of all complementarities
            all_pairs = [];
            all_products = [];
            std_indices = [];
            for k=1:problem_options.N_stages
                stage = obj.stages(k);
                for fe=stage.stage
                    all_pairs = [all_pairs;fe.all_comp_pairs];
                    std_indices = [std_indices;fe.ind_std_comp];
                end
            end
            obj.ind_std_comp = std_indices;
            obj.cross_comps = all_pairs;
            all_products = apply_psi(all_pairs, @(x,y,t) x.*y, 0);

            obj.comp_res = Function('comp_res', {obj.w, obj.p}, {max(all_products)});

            % Scalar-valued complementairity residual
            if problem_options.use_fesd
                J_comp_fesd = max(obj.cc_vector);
                J_comp = J_comp_fesd;
            else
                J_comp_fesd = J_comp_std;
                J_comp = J_comp_std;
            end

            % TODO Figure out if any of these are needed and cleanup how they are calculated
            %obj.comp_res = Function('comp_res', {obj.w, obj.p}, {J_comp});
            obj.comp_std = Function('comp_std', {obj.w, obj.p}, {J_comp_std});
            obj.comp_fesd = Function('comp_fesd', {obj.w, obj.p}, {J_comp_fesd});
            obj.augmented_objective_fun = Function('augmented_objective_fun', {obj.w, obj.p}, {obj.augmented_objective});
            obj.objective_fun = Function('objective_fun', {obj.w, obj.p}, {obj.objective});
            obj.g_fun = Function('g_fun', {obj.w, obj.p}, {obj.g});
            obj.G_fun = Function('G_fun', {obj.w}, {obj.cross_comps(:,1)});
            obj.H_fun = Function('H_fun', {obj.w}, {obj.cross_comps(:,2)});

            obj.p0 = [problem_options.rho_sot; problem_options.rho_h; problem_options.rho_terminal; problem_options.T];

            if dims.n_p_global > 0
                obj.p0 = [obj.p0; model.p_global_val];
            end

            if dims.n_p_time_var > 0
                obj.p0 = [obj.p0; model.p_time_var_val];
            end
            obj.w0_original = obj.w0;

            % Define CasADi function for the switch indicator function.
            nu_fun = Function('nu_fun', {obj.w,obj.p},{obj.nu_vector});
            obj.nu_fun = nu_fun;
            
            % create CasADi function for objective gradient.
            nabla_J = obj.augmented_objective.jacobian(obj.w);
            nabla_J_fun = Function('nabla_J_fun', {obj.w,obj.p},{nabla_J});
            obj.nabla_J = nabla_J;
            obj.nabla_J_fun = nabla_J_fun;

            mpcc_generating_time = toc;
            if problem_options.print_level >=2
                fprintf('MPCC generated in in %2.2f s. \n',mpcc_generating_time);
            end
        end

        % TODO this should be private
        function create_primal_variables(obj)
            import casadi.*
            fe0 = FiniteElementZero(obj.problem_options, obj.dims, obj.model);
            obj.fe0 = fe0;

            %             obj.p = vertcat(obj.p, fe0.x0, fe0.lambda{1,:},fe0.y_gap{1,:},fe0.gamma{1,:},fe0.gamma_d{1,:},fe0.delta_d{1,:},fe0.p_vt{1,:},fe0.n_vt{1,:});
            obj.p = vertcat(obj.p, fe0.x0, fe0.cross_comp_cont_0{1,:}, fe0.cross_comp_cont_1{1,:},fe0.cross_comp_cont_2{1,:});

            X0 = fe0.x{1};
            obj.addVariable(X0,...
                'x0',...
                fe0.lbw(fe0.ind_x{1}),...
                fe0.ubw(fe0.ind_x{1}),...
                fe0.w0(fe0.ind_x{1}));
            obj.addConstraint(fe0.g, fe0.lbg, fe0.ubg);
            prev_fe = fe0;

            s_sot = [];
            if obj.problem_options.time_rescaling && obj.problem_options.use_speed_of_time_variables
                if ~obj.problem_options.local_speed_of_time_variable
                    s_sot = define_casadi_symbolic(obj.problem_options.casadi_symbolic_mode, 's_sot', 1);
                    obj.addVariable(s_sot,...
                        'sot',...
                        obj.problem_options.s_sot_min,...
                        obj.problem_options.s_sot_max,...
                        obj.problem_options.s_sot0);
                    if obj.problem_options.time_freezing
                        obj.augmented_objective = obj.augmented_objective + obj.rho_sot_p*(s_sot-1)^2;
                    end
                end
            end

            for ii=1:obj.problem_options.N_stages
                % TODO: maybe this should be a function
                stage = ControlStage(prev_fe, obj.problem_options, obj.model, obj.dims, ii, s_sot, obj.T_final, obj.rho_h_p, obj.rho_sot_p);
                obj.stages = [obj.stages, stage];

                obj.addControlStage(stage);
                obj.augmented_objective = obj.augmented_objective + stage.augmented_objective;
                obj.objective = obj.objective + stage.objective;
                prev_fe = stage.stage(end);
            end
        end

        % TODO this should be private
        function addControlStage(obj, stage)
            w_len = length(obj.w);

            obj.addPrimalVector(stage.w, stage.lbw, stage.ubw, stage.w0);

            obj.ind_h = [obj.ind_h, increment_indices(stage.ind_h,w_len)];
            obj.ind_u = [obj.ind_u, {stage.ind_u+w_len}];
            obj.ind_sot = [obj.ind_sot, stage.ind_sot+w_len];
            obj.ind_x(stage.ctrl_idx, :, :) = increment_indices(stage.ind_x, w_len);
            obj.ind_v(stage.ctrl_idx, :, :) = increment_indices(stage.ind_v, w_len);
            obj.ind_theta(stage.ctrl_idx, :, :) = increment_indices(stage.ind_theta, w_len);
            obj.ind_lam(stage.ctrl_idx, :, :) = increment_indices(stage.ind_lam, w_len);
            obj.ind_mu(stage.ctrl_idx, :, :) = increment_indices(stage.ind_mu, w_len);
            obj.ind_alpha(stage.ctrl_idx, :, :) = increment_indices(stage.ind_alpha, w_len);
            obj.ind_lambda_n(stage.ctrl_idx, :, :) = increment_indices(stage.ind_lambda_n, w_len);
            obj.ind_lambda_p(stage.ctrl_idx, :, :) = increment_indices(stage.ind_lambda_p, w_len);
            obj.ind_beta(stage.ctrl_idx, :, :) = increment_indices(stage.ind_beta, w_len);
            obj.ind_theta_step(stage.ctrl_idx, :, :) = increment_indices(stage.ind_theta_step, w_len);
            obj.ind_z(stage.ctrl_idx, :, :) = increment_indices(stage.ind_z, w_len);
            obj.ind_nu_lift = [obj.ind_nu_lift, increment_indices(stage.ind_nu_lift, w_len)];

            obj.ind_lambda_normal(stage.ctrl_idx, :, :) = increment_indices(stage.ind_lambda_normal, w_len);
            obj.ind_lambda_tangent(stage.ctrl_idx, :, :) = increment_indices(stage.ind_lambda_tangent, w_len);
            obj.ind_y_gap(stage.ctrl_idx, :, :) = increment_indices(stage.ind_y_gap, w_len);
            obj.ind_gamma(stage.ctrl_idx, :, :) = increment_indices(stage.ind_gamma, w_len);
            obj.ind_beta_conic(stage.ctrl_idx, :, :) = increment_indices(stage.ind_beta_conic, w_len);
            obj.ind_gamma_d(stage.ctrl_idx, :, :) = increment_indices(stage.ind_gamma_d, w_len);
            obj.ind_beta_d(stage.ctrl_idx, :, :) = increment_indices(stage.ind_beta_d, w_len);
            obj.ind_delta_d(stage.ctrl_idx, :, :) = increment_indices(stage.ind_delta_d, w_len);
            obj.ind_p_vt(stage.ctrl_idx, :, :) = increment_indices(stage.ind_p_vt, w_len);
            obj.ind_n_vt(stage.ctrl_idx, :, :) = increment_indices(stage.ind_n_vt, w_len);
            obj.ind_alpha_vt(stage.ctrl_idx, :, :) = increment_indices(stage.ind_alpha_vt, w_len);

            obj.ind_x_left_bp(stage.ctrl_idx, :) = increment_indices(stage.ind_x_left_bp, w_len);
            obj.ind_Y_gap(stage.ctrl_idx, :) = increment_indices(stage.ind_Y_gap, w_len);
            obj.ind_Lambda_normal(stage.ctrl_idx, :) = increment_indices(stage.ind_Lambda_normal, w_len);
            obj.ind_Lambda_tangent(stage.ctrl_idx, :) = increment_indices(stage.ind_Lambda_tangent, w_len);
            obj.ind_Gamma(stage.ctrl_idx, :) = increment_indices(stage.ind_Gamma, w_len);
            obj.ind_Beta_conic(stage.ctrl_idx, :) = increment_indices(stage.ind_Beta_conic, w_len);
            obj.ind_Gamma_d(stage.ctrl_idx, :) = increment_indices(stage.ind_Gamma_d, w_len);
            obj.ind_Beta_d(stage.ctrl_idx, :) = increment_indices(stage.ind_Beta_d, w_len);
            obj.ind_Delta_d(stage.ctrl_idx, :) = increment_indices(stage.ind_Delta_d, w_len);
            obj.ind_L_vn(stage.ctrl_idx, :) = increment_indices(stage.ind_L_vn, w_len);
            obj.ind_P_vt(stage.ctrl_idx, :) = increment_indices(stage.ind_P_vt, w_len);
            obj.ind_N_vt(stage.ctrl_idx, :) = increment_indices(stage.ind_N_vt, w_len);
            obj.ind_Alpha_vt(stage.ctrl_idx, :) = increment_indices(stage.ind_Alpha_vt, w_len);

            obj.ind_comp_lift = [obj.ind_comp_lift; increment_indices(stage.ind_comp_lift,w_len)];

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

        function cc_vector = get.cc_vector(obj)
            cc_vector = [];
            for stage=obj.stages
                for fe=stage.stage
                    cc_vector = vertcat(cc_vector, fe.all_comp_pairs .* fe.all_comp_pairs);
                end
            end
        end

        function u = get.u(obj)
            u = cellfun(@(u) obj.w(u), obj.ind_u, 'UniformOutput', false);
        end

        function sot = get.sot(obj)
            sot = cellfun(@(sot) obj.w(sot), obj.ind_sot, 'UniformOutput', false);
        end

        function nu_vector = get.nu_vector(obj)
            nu_vector = [];
            for k=obj.problem_options.N_stages
                stage = obj.stages(k);
                for fe=stage.stage
                    nu_vector = vertcat(nu_vector,fe.nu_vector);
                end
            end
        end

        function ind_z_all = get.ind_z_all(obj)
            ind_z_all = [flatten_ind(obj.ind_theta(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_lam(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_mu(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_alpha(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_lambda_n(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_lambda_p(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_beta(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_theta_step(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_z(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_lambda_normal(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_lambda_tangent(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_y_gap(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_gamma(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_beta_conic(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_gamma_d(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_beta_d(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_delta_d(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_p_vt(:,:,1:obj.dims.n_s))
                flatten_ind(obj.ind_p_vt(:,:,1:obj.dims.n_s))
                        ];
            ind_z_all = sort(ind_z_all);
        end

        function ind_x_all = get.ind_x_all(obj)
            ind_x_all = [obj.ind_x0.'; flatten_ind(obj.ind_x)];
        end

        function print(obj,filename)
            if exist('filename')
                delete(filename);
                fileID = fopen(filename, 'w');
            else
                fileID = 1;
            end
            fprintf(fileID, "i\tlbg\t\t ubg\t\t g_expr\n");
            for ii = 1:length(obj.lbg)
                expr_str = formattedDisplayText(obj.g(ii));
                fprintf(fileID, "%d\t%.2e\t%.2e\t%s\n", ii, obj.lbg(ii), obj.ubg(ii), expr_str);
            end

            fprintf(fileID, "\nw\t\t\tw0\t\tlbw\t\tubw\n");
            for ii = 1:length(obj.lbw)
                % keyboard
                expr_str = pad(formattedDisplayText(obj.w(ii)), 20);
                lb_str = pad(sprintf('%.2e', obj.lbw(ii)), 10);
                fprintf(fileID, "%s\t%.2e\t%s\t%.2e\t\n", expr_str, obj.w0(ii), lb_str, obj.ubw(ii));
            end

            fprintf(fileID, "\nCross Complementarity Pairs\n");
            fprintf(fileID, "\na \t\t\t\t\t\t b\n");

            for stage=obj.stages
                for fe=stage.stage
                    for ii=1:size(fe.all_comp_pairs, 1)
                        a_str = pad(formattedDisplayText(fe.all_comp_pairs(ii,1)), 20);
                        b_str = pad(formattedDisplayText(fe.all_comp_pairs(ii,2)), 20);
                        fprintf(fileID, "%s\t\t\t\t%s\n", a_str, b_str);
                    end
                end
            end

            fprintf(fileID, '\naugmented objective\n');
            fprintf(fileID, strcat(formattedDisplayText(obj.augmented_objective), '\n'));
        end

        function json = jsonencode(obj,varargin)
            import casadi.*
            mpcc_struct = struct(obj);

            mpcc_struct.model = obj.model;
            mpcc_struct.problem_options = obj.problem_options

            json = jsonencode(mpcc_struct);
        end
    end

    methods(Static)
        function obj = from_json(json)
            mpcc_struct = jsondecode(json);
            model = NosnocModel.from_struct(mpcc_struct.model);
            problem_options = NosnocProblemOptions.from_struct(mpcc_struct.problem_options);
            obj = NosnocMpcc(problem_options, model.dims, model);
        end
    end
end
