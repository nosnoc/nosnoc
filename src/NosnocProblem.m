% BSD 2-Clause License

% Copyright (c) 2023, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

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


% TODO cleanup steps:
%      - Create primal variables all at once.
%      - Separate sections into separate functions operating on the `problem` struct/class
%      - time variables should probably not just be lumped into the state, for readability.
%      - remove index in symbolic variable defintions and add instructive
%        names, e.g., Uk -> U,  h_ki -> h_fe, X_ki_stages ->  X_rk_stages
%      - provide instructive names for terminal constraint relaxations
%      - provide more instructive names for cross_comp (match python)

classdef NosnocProblem < NosnocFormulationObject
    properties
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
%         ind_P_vn
%         ind_N_vn
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

        % Parameter index variables
        ind_p_x0
        ind_p_global
        ind_p_time_var

        % Problem data
        model
        settings
        dims
        ocp

        % original initialization
        w0_original

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
        fe0 % Zeroth finite element (contains X0, lambda00)
        stages % control stages

        % complementarity residual functions
        comp_res
        comp_std
        comp_fesd

        % Problem cost function
        cost_fun

        % Problem objective function
        objective_fun

        % Problem constraint function
        g_fun

        % switch indicator function
        nu_fun

        nabla_J
        nabla_J_fun
    end
    % remaining list of TODOs
    % TODO: cleanup/add properties (in all components)
    % TODO: Create solver object, which will interact with setting parameters.

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
        function obj = NosnocProblem(settings, dims, model)
            import casadi.*
            obj@NosnocFormulationObject();

            if settings.right_boundary_point_explicit || settings.dcs_mode == DcsMode.CLS
                rbp_allowance = 0;
            else
                rbp_allowance = 1;
            end

            if settings.dcs_mode == "CLS" && ~settings.right_boundary_point_explicit
                rbp_x_only = 1;
            else
                rbp_x_only = 0;
            end

            right_ygap = 0;
            if settings.dcs_mode == "CLS" && ~settings.right_boundary_point_explicit
                right_ygap = 1;
            end

            obj.ind_u = [];
            obj.ind_x0 = [];
            obj.ind_x = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance+rbp_x_only);
            obj.ind_v = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s);
            obj.ind_z = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            % Stewart
            obj.ind_theta = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_lam = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_mu = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            % Step
            obj.ind_alpha = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_lambda_n = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_lambda_p = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_beta = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_theta_step = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            % CLS
            obj.ind_lambda_normal = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_lambda_tangent = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_y_gap = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+right_ygap+rbp_allowance);
            % friction multipliers and lifting
            % conic
            obj.ind_gamma =  cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_beta_conic = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            % poly
            obj.ind_gamma_d = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_beta_d = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_delta_d = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            % variables related to conic
            obj.ind_p_vt = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_n_vt = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);
            obj.ind_alpha_vt = cell(settings.N_stages,settings.N_finite_elements(1),dims.n_s+rbp_allowance);

            % Impulse
            obj.ind_x_left_bp = cell(settings.N_stages,settings.N_finite_elements(1),1);
            obj.ind_Y_gap = cell(settings.N_stages,settings.N_finite_elements(1),1);
            obj.ind_Lambda_normal = cell(settings.N_stages,settings.N_finite_elements(1),1);
            obj.ind_Lambda_tangent = cell(settings.N_stages,settings.N_finite_elements(1),1);
            obj.ind_Gamma = cell(settings.N_stages,settings.N_finite_elements(1),1);
            obj.ind_Beta_conic = cell(settings.N_stages,settings.N_finite_elements(1),1);
            obj.ind_Gamma_d = cell(settings.N_stages,settings.N_finite_elements(1),1);
            obj.ind_Beta_d = cell(settings.N_stages,settings.N_finite_elements(1),1);
            obj.ind_Delta_d = cell(settings.N_stages,settings.N_finite_elements(1),1);
%             obj.ind_P_vn = cell(settings.N_stages,settings.N_finite_elements(1),1);
%             obj.ind_N_vn = cell(settings.N_stages,settings.N_finite_elements(1),1);
            obj.ind_L_vn = cell(settings.N_stages,settings.N_finite_elements(1),1);
            obj.ind_P_vt = cell(settings.N_stages,settings.N_finite_elements(1),1);
            obj.ind_N_vt = cell(settings.N_stages,settings.N_finite_elements(1),1);
            obj.ind_Alpha_vt = cell(settings.N_stages,settings.N_finite_elements(1),1);

            
            % misc
            obj.ind_nu_lift = {};
            obj.ind_h = {};
            obj.ind_sot = {};
            obj.ind_v_global = [];

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
                obj.ind_p_global = n_p+1:n_p+dims.n_p_global;
                obj.p = [obj.p; model.p_global];
            end
            if dims.n_p_time_var > 0;
                n_p = length(obj.p);
                obj.ind_p_time_var = arrayfun(@(s) (n_p+(s*dims.n_p_time_var)+1):(n_p+(s*dims.n_p_time_var)+dims.n_p_time_var) , 0:settings.N_stages-1);
                obj.p = [obj.p; model.p_time_var_stages(:)];
            end
            if settings.time_optimal_problem
                % the final time in time optimal control problems
                T_final = define_casadi_symbolic(settings.casadi_symbolic_mode, 'T_final', 1);
                obj.T_final = T_final;
                T_final_guess = model.T;
            end

            % Add global vars
            obj.addVariable(model.v_global,...
                'v_global',...
                model.lbv_global,...
                model.ubv_global,...
                model.v0_global)

            obj.create_primal_variables();

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
                            for k=1:settings.N_stages
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
                        for k=1:settings.N_stages
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
                                for k=1:settings.N_stages
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
            if ~isempty(model.g_terminal)
                if settings.relax_terminal_constraint_homotopy
                    rho_terminal_p = 1/sigma_p;
                end
                X_end = last_fe.x{end};
                g_terminal = model.g_terminal_fun(X_end, model.p_global, model.v_global);
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
            obj.cost = obj.cost + model.f_lsq_T_fun(last_fe.x{end}, model.x_ref_end_val, model.p_global);
            obj.objective = obj.objective + model.f_lsq_T_fun(last_fe.x{end}, model.x_ref_end_val, model.p_global);

            % Process terminal costs
            obj.cost = obj.cost + model.f_q_T_fun(last_fe.x{end}, model.p_global, model.v_global);
            obj.objective = obj.objective + model.f_q_T_fun(last_fe.x{end}, model.p_global, model.v_global);

            % Process elastic costs
            if settings.elasticity_mode == ElasticityMode.ELL_INF
                if settings.objective_scaling_direct
                    obj.cost = obj.cost + (1/sigma_p)*obj.s_elastic;
                else
                    obj.cost = sigma_p*obj.cost + obj.s_elastic;
                end
            end
            if settings.elasticity_mode == ElasticityMode.ELL_1
                sum_s_elastic = 0;
                for k=1:settings.N_stages
                    stage=obj.stages(k);
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
                obj.objective = obj.objective + T_final;
            end

            % Calculate standard complementarities.
            J_comp_std = 0;
            J_comp_std_infty = 0;
            for k=1:settings.N_stages
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
            for k=1:settings.N_stages
                stage = obj.stages(k);
                for fe=stage.stage
                    all_pairs = [all_pairs;fe.all_comp_pairs];
                end
            end
            all_products = apply_psi(all_pairs, @(x,y,t) x*y, 0);

            obj.comp_res = Function('comp_res', {obj.w, obj.p}, {max(all_products)});

            % Scalar-valued complementairity residual
            if settings.use_fesd
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
            obj.cost_fun = Function('cost_fun', {obj.w, obj.p}, {obj.cost});
            obj.objective_fun = Function('objective_fun', {obj.w, obj.p}, {obj.objective});
            obj.g_fun = Function('g_fun', {obj.w, obj.p}, {obj.g});

            obj.p0 = [settings.sigma_0; settings.rho_sot; settings.rho_h; settings.rho_terminal; model.T];

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
            nabla_J = obj.cost.jacobian(obj.w);
            nabla_J_fun = Function('nabla_J_fun', {obj.w,obj.p},{nabla_J});
            obj.nabla_J = nabla_J;
            obj.nabla_J_fun = nabla_J_fun;
        end

        % TODO this should be private
        function create_primal_variables(obj)
            import casadi.*
            fe0 = FiniteElementZero(obj.settings, obj.dims, obj.model);
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

            if obj.settings.elasticity_mode == ElasticityMode.ELL_INF
                s_elastic = define_casadi_symbolic(obj.settings.casadi_symbolic_mode, 's_elastic',1);
                obj.s_elastic = s_elastic;
                if obj.settings.elastic_scholtes
                    obj.settings.s_elastic_max = inf;
                    obj.addConstraint(s_elastic-obj.sigma_p,-inf,0);
                end
                obj.addVariable(s_elastic, 'elastic', obj.settings.s_elastic_min, obj.settings.s_elastic_max, obj.settings.s_elastic_0);
            else
                s_elastic = [];
            end

            for ii=1:obj.settings.N_stages
                % TODO: maybe this should be a function
                stage = ControlStage(prev_fe, obj.settings, obj.model, obj.dims, ii, s_sot, obj.T_final, obj.sigma_p, obj.rho_h_p, obj.rho_sot_p, s_elastic);
                obj.stages = [obj.stages, stage];

                obj.addControlStage(stage);
                obj.cost = obj.cost + stage.cost;
                obj.objective = obj.objective + stage.objective;
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

            psi_fun = settings.psi_fun;

            if settings.elasticity_mode == ElasticityMode.NONE
                sigma = sigma_p;
            else
                sigma = s_elastic;
            end

            g_cross_comp = SX([]);
            % TODO Implement other modes
            if ~settings.use_fesd || settings.cross_comp_mode < 11
                % Do nothing, handled at the FE or stage level
                return
            elseif settings.cross_comp_mode == 11
                for r=1:dims.n_sys
                    g_r = 0;
                    nz_r = [];
                    for stage=obj.stages
                        for fe=stage.stage
                            pairs = fe.cross_comp_pairs(:, :, r);
                            expr_cell = cellfun(@(pair) apply_psi(pair, psi_fun, sigma), pairs, 'uni', false);
                            nonzeros = cellfun(@(x) vector_is_zero(x), expr_cell, 'uni', 0);
                            if size(vertcat(expr_cell{:}), 1) == 0
                                exprs= [];
                            elseif settings.relaxation_method == RelaxationMode.TWO_SIDED
                                exprs_p = cellfun(@(c) c(:,1), expr_cell, 'uni', false);
                                exprs_n = cellfun(@(c) c(:,2), expr_cell, 'uni', false);
                                nonzeros_p = cellfun(@(x) vector_is_zero(x), exprs_p, 'uni', 0);
                                nonzeros_n = cellfun(@(x) vector_is_zero(x), exprs_n, 'uni', 0);
                                nonzeros = [sum([nonzeros_p{:}], 2),sum([nonzeros_n{:}], 2)]';
                                exprs = [sum2([exprs_p{:}]),sum2([exprs_n{:}])]';
                                exprs = exprs(:);
                            else
                                nonzeros = sum([nonzeros{:}], 2);
                                exprs = sum2([expr_cell{:}]);
                            end
                            if isempty(nz_r)
                                nz_r = zeros(size(nonzeros));
                            end
                            g_r = g_r + extract_nonzeros_from_vector(exprs);
                            nz_r = nz_r + nonzeros;
                        end
                    end
                    g_r = scale_sigma(g_r, sigma, nz_r);
                    g_cross_comp = vertcat(g_cross_comp, g_r);
                end
            elseif settings.cross_comp_mode == 12
                for r=1:dims.n_sys
                    g_r = 0;
                    nz_r = [];
                    for stage=obj.stages
                        for fe=stage.stage
                            pairs = fe.cross_comp_pairs(:, :, r);
                            expr_cell = cellfun(@(pair) apply_psi(pair, psi_fun, sigma), pairs, 'uni', false);
                            nonzeros = cellfun(@(x) vector_is_zero(x), expr_cell, 'uni', 0);
                            if size(vertcat(expr_cell{:}), 1) == 0
                                exprs= [];
                            elseif settings.relaxation_method == RelaxationMode.TWO_SIDED
                                exprs_p = cellfun(@(c) c(:,1), expr_cell, 'uni', false);
                                exprs_n = cellfun(@(c) c(:,2), expr_cell, 'uni', false);
                                nonzeros_p = cellfun(@(x) vector_is_zero(x), exprs_p, 'uni', 0);
                                nonzeros_n = cellfun(@(x) vector_is_zero(x), exprs_n, 'uni', 0);
                                nonzeros = [sum(sum([nonzeros_p{:}], 2), 1),sum(sum([nonzeros_n{:}], 2), 1)]';
                                exprs = [sum1(sum2([exprs_p{:}])),sum1(sum2([exprs_n{:}]))]';
                                exprs = exprs(:);
                            else
                                nonzeros = sum(sum([nonzeros{:}], 2),1);
                                exprs = sum1(sum2([expr_cell{:}]));
                            end
                            if isempty(nz_r)
                                nz_r = zeros(size(nonzeros));
                            end
                            g_r = g_r + extract_nonzeros_from_vector(exprs);
                            nz_r = nz_r + nonzeros;
                        end
                    end
                    g_r = scale_sigma(g_r, sigma, nz_r);
                    g_cross_comp = vertcat(g_cross_comp, g_r);
                end
            end

            % If We need to add a cost from the reformulation do that as needed;
            if settings.mpcc_mode == MpccMode.ell_1_penalty % this implies bilinear
                cost = 0;
                for r=1:dims.n_sys
                    for stage=obj.stages
                        for fe=stage.stage
                            pairs = fe.cross_comp_pairs(:, :, r);
                            expr_cell = cellfun(@(pair) apply_psi(pair, @(a,b,t) a*b, 0), pairs, 'uni', false);
                            expr = sum1(sum2([expr_cell{:}]));
                            cost = cost + expr;
                        end
                    end
                end
                if settings.objective_scaling_direct
                    obj.cost = obj.cost + (1/sigma_p)*cost;
                else
                    obj.cost = sigma_p*obj.cost + cost;
                end
            else
                g_comp = g_cross_comp;
                n_comp = length(g_cross_comp);

                [g_comp_lb, g_comp_ub, g_comp] = generate_mpcc_relaxation_bounds(g_comp, settings);

                % Add reformulated constraints
                obj.addConstraint(g_comp, g_comp_lb, g_comp_ub);

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
%             obj.ind_P_vn(stage.ctrl_idx, :) = increment_indices(stage.ind_P_vn, w_len);
%             obj.ind_N_vn(stage.ctrl_idx, :) = increment_indices(stage.ind_N_vn, w_len);
            obj.ind_L_vn(stage.ctrl_idx, :) = increment_indices(stage.ind_L_vn, w_len);
            obj.ind_P_vt(stage.ctrl_idx, :) = increment_indices(stage.ind_P_vt, w_len);
            obj.ind_N_vt(stage.ctrl_idx, :) = increment_indices(stage.ind_N_vt, w_len);
            obj.ind_Alpha_vt(stage.ctrl_idx, :) = increment_indices(stage.ind_Alpha_vt, w_len);

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
                    for ic = 1:size(fe.all_comp_pairs, 1)
                        cc_vector = vertcat(cc_vector, fe.all_comp_pairs(ic, 1) * fe.all_comp_pairs(ic, 2));
                    end
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
            for k=obj.settings.N_stages
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
            for i = 1:length(obj.lbg)
                expr_str = formattedDisplayText(obj.g(i));
                fprintf(fileID, "%d\t%.2e\t%.2e\t%s\n", i, obj.lbg(i), obj.ubg(i), expr_str);
            end

            fprintf(fileID, "\nw\t\t\tw0\t\tlbw\t\tubw\n");
            for i = 1:length(obj.lbw)
                % keyboard
                expr_str = pad(formattedDisplayText(obj.w(i)), 20);
                lb_str = pad(sprintf('%.2e', obj.lbw(i)), 10);
                fprintf(fileID, "%s\t%.2e\t%s\t%.2e\t\n", expr_str, obj.w0(i), lb_str, obj.ubw(i));
            end

            fprintf(fileID, '\nobjective\n');
            fprintf(fileID, strcat(formattedDisplayText(obj.cost), '\n'));
        end
    end
end
