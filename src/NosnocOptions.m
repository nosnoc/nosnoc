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
classdef NosnocOptions < handle
% TODO clean up much of the work here.
    properties
        % General
        solver_name {mustBeTextScalar} = 'nosnoc_solver'
        nlpsol = 'ipopt'
        use_fesd(1,1) logical = 1
        casadi_symbolic_mode {mustBeMember(casadi_symbolic_mode,{'casadi.SX', 'casadi.MX'})} = 'casadi.SX'

        % descritization
        N_stages(1,1) {mustBeInteger, mustBePositive} = 1;
        N_finite_elements {mustBeInteger, mustBePositive} = 2;

        % Integrator settings
        real_time_plot(1,1) logical = 0
        break_simulation_if_infeasible(1,1) logical = 1

        % IRK and FESD Settings
        n_s(1,1) {mustBeInteger, mustBeInRange(n_s, 1, 9)} = 2
        irk_scheme(1,1) IRKSchemes = IRKSchemes.RADAU_IIA
        irk_representation IrkRepresentation = IrkRepresentation.integral;
        cross_comp_mode {mustBeInRange(cross_comp_mode, 1, 12)}= 3
        gamma_h(1,1) double {mustBeReal, mustBeInRange(gamma_h, 0, 1)} = 1
        dcs_mode DcsMode = DcsMode.Stewart

        % Initialization - Stewart
        lp_initialization(1,1) logical = 0
        initial_theta(1,1) double {mustBeReal, mustBeFinite} = 1
        initial_lambda(1,1) double {mustBeReal, mustBeFinite} = 1
        initial_mu(1,1) double {mustBeReal, mustBeFinite} = 1

        % Initialization - Step
        initial_alpha(1,1) double {mustBeReal, mustBeFinite} = 1
        initial_lambda_0(1,1) double {mustBeReal, mustBeFinite} = 1
        initial_lambda_1(1,1) double {mustBeReal, mustBeFinite} = 1
        initial_beta(1,1) double {mustBeReal, mustBeFinite} = 1
        initial_theta_step(1,1) double {mustBeReal, mustBeFinite} = 1

        % Step theta lifting
        pss_lift_step_functions(1,1) logical = 0
        n_depth_step_lifting(1,1) {mustBeInteger, mustBeGreaterThanOrEqual(n_depth_step_lifting, 2)} = 2

        %General NLP/OCP Settings
        g_path_at_fe(1,1) logical = 0 % evaluate nonlinear path constraint at every finte element boundary
        g_path_at_stg(1,1) logical = 0 % evaluate nonlinear path constraint at every stage

        x_box_at_fe(1,1) logical = 1 % evaluate box constraint for diff states at every finite element boundary point
        x_box_at_stg(1,1) logical = 1 % evaluate box constraint for diff states at every stage point. (is set to zero per default in differential irk mode, as it becomes a linear instead of box constraint)

        time_optimal_problem(1,1) = 0
        simple_v0_guess(1,1) = 0 % TODO what is this

        % MPCC and Homotopy Settings
        comp_tol(1,1) double {mustBeReal, mustBePositive} = 1e-9
        mpcc_mode(1,1) MpccMode = MpccMode.Scholtes_ineq % 'direct', 'Scholtes_eq', 'Scholtes_ineq', 'ell_1_penalty', 'elastic_ineq', 'elastic_eq' , 'elastic_two_sided',
                                           % 'elastic_ell_1_ineq', 'elastic_ell_1_eq', 'elastic_ell_1_two_sided'
        objective_scaling_direct(1,1) logical = 1
        sigma_0(1,1) double {mustBeReal, mustBePositive} = 1
        sigma_N(1,1) double {mustBeReal, mustBePositive} = 1e-9
        homotopy_update_rule = 'linear' % 'linear' sigma_k = homotopy_update_slope*sigma_N
                                        % 'superlinear' - sigma_k = max(sigma_N,min(homotopy_update_slope*sigma_k,sigma_k^homotopy_update_exponent))
                                        % TODO enum
        homotopy_update_slope(1,1) double {mustBeReal, mustBeInRange(homotopy_update_slope, 0, 1, 'exclusive')} = 0.1
        homotopy_update_exponent(1,1) double {mustBeReal, mustBePositive} = 1.5 % the exponent in the superlinear rule
        N_homotopy = 0 % 0 -> set automatically
        s_elastic_max(1,1) double {mustBeReal, mustBePositive} = 1e1
        s_elastic_min(1,1) double {mustBeReal, mustBeNonnegative} = 0
        s_elastic_0(1,1) double {mustBeReal, mustBePositive} = 1
        elastic_scholtes(1,1) logical = 0

        % Default settings for the barrier tuned penalty/slack variables for mpcc modes 8 do 10.
        rho_penalty(1,1) double {mustBeReal, mustBePositive} = 1e1
        sigma_penalty = 0 % TODO what is this

        rho_lambda(1,1) double {mustBeReal, mustBePositive} = 1
        rho_scale(1,1) double {mustBeReal, mustBePositive} = 30

        sigma_scale(1,1) double {mustBeReal, mustBePositive} = 0.1

        rho_min(1,1) double {mustBeReal, mustBePositive} = 0.1
        rho_max(1,1) double {mustBeReal, mustBePositive} = (log(30)-log(1e-16))/1

        rho_0(1,1) double {mustBeReal, mustBePositive} = 0.5
        % sigma_0 = sigma_scale*rho_scale*exp(-rho_lambda*rho_0)
        nonlinear_sigma_rho_constraint(1,1) logical = 1
        convex_sigma_rho_constraint(1,1) logical = 0

        ratio_for_homotopy_stop(1,1) double {mustBeReal, mustBePositive} = 0.75

        % Homotopy preprocess and polishing steps
%         h_fixed_iterations(1,1) logical = 0
%         h_fixed_max_iter(1,1) logical = 1 % number of iterations that are done with fixed h in the homotopy loop
%         h_fixed_change_sigma(1,1) logical = 1 % if this is on, do not update sigma and just solve on nlp with fixed h.
        polishing_step(1,1) logical = 0 % heuristic for fixing active set, yet exerimental, not recommended to use.
        polishing_derivative_test(1,1) logical = 0 % check in sliding mode also the derivative of switching functions
        h_fixed_to_free_homotopy(1,1) logical = 0 % start with large penaly for equidistant grid, end with variable equilibrated grid.


        % Step equilibration
        rho_h(1,1) double {mustBeReal, mustBePositive} = 1
        step_equilibration(1,1) StepEquilibrationMode = StepEquilibrationMode.heuristic_mean % heuristic_mean, l2_relaxed, l2_relaxed_scaled, direct, direct_homotopy, direct_homotopy_lift
        step_equilibration_sigma(1,1) double {mustBeReal, mustBePositive} = 0.1 % slope at zero in rescaling the indicator function, nu_ki_rescaled = tanh(nu_ki/step_equilibration_sigma)

        % Multiple shooting type discretization
        equidistant_control_grid(1,1) logical = 1

        % Time-Setting & Time-Freezing
        time_freezing(1,1) logical = 0
        time_freezing_inelastic(1,1) logical = 0
        % for time optimal problems and equidistant control grids in physical time
        use_speed_of_time_variables(1,1) logical = 0
        local_speed_of_time_variable(1,1) logical = 0
        stagewise_clock_constraint(1,1) logical = 1
        impose_terminal_phyisical_time(1,1) logical = 1
        s_sot0(1,1) double {mustBeReal, mustBePositive} = 1
        s_sot_max(1,1) double {mustBeReal, mustBePositive} = 25
        s_sot_min(1,1) double {mustBeReal, mustBePositive} = 1
        S_sot_nominal(1,1) double {mustBeReal, mustBePositive} = 1
        rho_sot(1,1) double {mustBeReal, mustBeNonnegative} = 0

        T_final_max(1,1) double {mustBeReal, mustBePositive} = 1e2
        T_final_min(1,1) double {mustBeReal, mustBeNonnegative} = 0
        time_freezing_reduced_model(1,1) logical = 0 % analytic reduction of lifter formulation, less algebraic variables (experimental)
        time_freezing_hysteresis(1,1) logical = 0 % do not do automatic time freezing generation for hysteresis, it is not supported yet.
        time_freezing_nonlinear_friction_cone(1,1) logical = 1 % 1 - use nonlienar friction cone, 0 - use polyhedral l_inf approximation.

        time_freezing_quadrature_state(1,1) logical = 0 % make a nonsmooth quadrature state to integrate only if physical time is running
        time_freezing_lift_forces(1,1) logical = 0 % replace \dot{v} = M(q)^{-1}f(q,v,u) by dot{v} = z,  M(q)z - f(q,v,u) = 0

        % exerimental expert options
        nonsmooth_switching_fun(1,1) logical = 0 % experimental: use c = max(c1,c2) insetad of c = [c1c2]
        % expert mode: stabilize auxiliary dynamics in \nabla f_c(q) direction
        stabilizing_q_dynamics(1,1) logical = 0
        kappa_stabilizing_q_dynamics(1,1) double {mustBeReal, mustBePositive} = 1e-5
        % Verbose
        print_level = 3
        print_details_if_infeasible = 0;
        pause_homotopy_solver_if_infeasible = 0;

        % Settings specific to CLS
        friction_model (1,1) FrictionModel = FrictionModel.Conic; % use nonlinear friction ('Conic') or polyhedral approximation ('Polyhedral');
        conic_model_switch_handling (1,1) ConicModelSwitchHandling = ConicModelSwitchHandling.Abs; % How to treat switch detection in v_t.
        kappa_friction_reg (1,1) double {mustBeReal, mustBeNonnegative} = 0; % reg. term in friction equations to avoid large multipliers if no contact happens.
        lift_velocity_state(1,1) logical = 0; % define auxliary algebraic vairable, dot{v} = z_v, to avoid symbolic inversion of the inertia matrix;
        eps_cls double = 1e-3 % constraint: enforce f_c at Euler step with h * eps_cls

        % IPOPT Settings
        tol_ipopt(1,1) double {mustBeReal, mustBePositive} = 1e-16
        opts_casadi_nlp

        % Relaxation of terminal constraint
        relax_terminal_constraint(1,1) {mustBeInteger, mustBeInRange(relax_terminal_constraint, 0, 3)} = 0 %  0  - hard constraint, 1 - ell_1 , 2  - ell_2 , 3 - ell_inf TODO enum
        relax_terminal_constraint_from_above(1,1) logical = 0
        rho_terminal(1,1) double {mustBeReal, mustBePositive} = 1e2
        relax_terminal_constraint_homotopy(1,1) logical = 0 % terminal penalty is governed by homotopy parameter

        % Integrator Specific
        use_previous_solution_as_initial_guess(1,1) logical = 0
        simulation_problem(1,1) logical = 0

        % Misc TODO these should probably live in model
        general_inclusion(1,1) logical = 0
        there_exist_free_x0(1,1) logical = 0
        time_freezing_model_exists(1,1) logical = 0

        % TODO: make proper multiple solver class.
        multiple_solvers(1,1) logical = 0

        % Experimental:
        no_initial_impacts(1,1) logical = 0

        % All NLP parameters
        T_val(1,1) double {mustBeReal, mustBePositive} = 1
        p_val

        % Butcher Tableu
        A_irk double
        B_irk double
        b_irk double
        C_irk double
        D_irk double
        c_irk double

        % psi func
        psi_fun_type CFunctionType = CFunctionType.BILINEAR
        relaxation_method(1,1) RelaxationMode = RelaxationMode.INEQ
        elasticity_mode(1,1) ElasticityMode = ElasticityMode.NONE
        psi_fun

        right_boundary_point_explicit(1,1) logical % TODO this shoud live in model probably

        % Output options
        store_integrator_step_results(1,1) logical = 0
        ipopt_callback = [] % This should be a function handle that takes (model,problem,settings,ipopt_solver,results)
    end

    properties(Dependent)
        time_rescaling
        
    end

    methods
        function obj = NosnocOptions()

            default_tol = 1e-12;

            obj.opts_casadi_nlp.ipopt.print_level = 0;
            obj.opts_casadi_nlp.print_time = 0;
            obj.opts_casadi_nlp.ipopt.sb = 'yes';
            obj.opts_casadi_nlp.verbose = false;
            obj.opts_casadi_nlp.ipopt.max_iter = 500;
            obj.opts_casadi_nlp.ipopt.bound_relax_factor = 0;
            obj.opts_casadi_nlp.ipopt.tol = default_tol;
            obj.opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
            obj.opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
            obj.opts_casadi_nlp.ipopt.compl_inf_tol = default_tol;
            obj.opts_casadi_nlp.ipopt.mu_strategy = 'adaptive';
            obj.opts_casadi_nlp.ipopt.mu_oracle = 'quality-function';
            obj.opts_casadi_nlp.snopt = struct();

            obj.p_val = [obj.sigma_0,obj.rho_sot,obj.rho_h,obj.rho_terminal,obj.T_val];
        end

        function [] = preprocess(obj)
            import casadi.*
            % automatically set up casadi opts
            if obj.print_level < 4 
                obj.opts_casadi_nlp.ipopt.print_level=0;
                obj.opts_casadi_nlp.print_time=0;
                obj.opts_casadi_nlp.ipopt.sb= 'yes';
            elseif print_level == 4
                obj.opts_casadi_nlp.ipopt.print_level=0;
                obj.opts_casadi_nlp.print_time=1;
                obj.opts_casadi_nlp.ipopt.sb= 'no';
            else
                obj.opts_casadi_nlp.ipopt.print_level = 5;
            end

            % check irk scheme compatibility
            if ismember(obj.irk_scheme, IRKSchemes.differential_only)
                if obj.print_level >=1
                    fprintf(['Info: The user provided RK scheme: ' char(obj.irk_scheme) ' is only available in the differential representation.\n']);
                end
                obj.irk_representation = 'differential';
            end

            % check impacts mode
            if obj.dcs_mode == DcsMode.CLS && obj.time_freezing == 1
                if obj.print_level >=1
                    fprintf(['nosnoc: User uses dcs_mode.CLS and time_freezing = true at the same time. Defaulting to time freezing as it is currently more competitive\n']);
                end
                obj.time_freezing = 1;
                obj.dcs_mode = DcsMode.Step;
            end

            if obj.time_freezing
                obj.use_speed_of_time_variables = 1;
                obj.local_speed_of_time_variable = 1;
            end

            % N finite elements shape 
            if isscalar(obj.N_finite_elements)
                obj.N_finite_elements = obj.N_finite_elements*ones(obj.N_stages,1);
            else
                if length(obj.N_finite_elements) == obj.N_stages
                    obj.N_finite_elements = obj.N_finite_elements;
                else
                    error('settings.N_finite_elements must be length 1 or N_stages');
                end
            end

            % MPCC mode setup
            if obj.mpcc_mode == MpccMode.direct
                % TODO maybe need to check later if people try to set sigmas/n_homotopy as for now we just let them break the
                %      assumptions.
                if obj.print_level >= 1
                    fprintf('Info: Setting N_homotopy to 1 and sigma_0,sigma_N to constants in direct mode.\n')
                end
                obj.N_homotopy = 1;
                obj.sigma_0 = 1e-12;
                obj.sigma_N = 1e-12;
                obj.mpcc_mode = MpccMode.Scholtes_ineq;
            end
            if obj.mpcc_mode == MpccMode.ell_1_penalty
                if obj.print_level >= 1
                    fprintf('Info: Setting cross_comp_mode to 12 in ell_1_penalty mode.\n')
                end
                obj.cross_comp_mode = 12;
            end

            switch obj.mpcc_mode
                case MpccMode.Scholtes_ineq
                    obj.psi_fun_type = CFunctionType.BILINEAR;
                    obj.relaxation_method = RelaxationMode.INEQ;
                    obj.elasticity_mode = ElasticityMode.NONE;
                case MpccMode.Scholtes_eq
                    obj.psi_fun_type = CFunctionType.BILINEAR;
                    obj.relaxation_method = RelaxationMode.EQ;
                    obj.elasticity_mode = ElasticityMode.NONE;
                case MpccMode.elastic_ineq
                    obj.psi_fun_type = CFunctionType.BILINEAR;
                    obj.relaxation_method = RelaxationMode.INEQ;
                    obj.elasticity_mode = ElasticityMode.ELL_INF;
                case MpccMode.elastic_eq
                    obj.psi_fun_type = CFunctionType.BILINEAR;
                    obj.relaxation_method = RelaxationMode.EQ;
                    obj.elasticity_mode = ElasticityMode.ELL_INF;
                case MpccMode.elastic_two_sided
                    obj.psi_fun_type = CFunctionType.BILINEAR_TWO_SIDED;
                    obj.relaxation_method = RelaxationMode.TWO_SIDED;
                    obj.elasticity_mode = ElasticityMode.ELL_INF;
                case MpccMode.elastic_ell_1_ineq
                    obj.psi_fun_type = CFunctionType.BILINEAR;
                    obj.relaxation_method = RelaxationMode.INEQ;
                    obj.elasticity_mode = ElasticityMode.ELL_1;
                case MpccMode.elastic_ell_1_eq
                    obj.psi_fun_type = CFunctionType.BILINEAR;
                    obj.relaxation_method = RelaxationMode.EQ;
                    obj.elasticity_mode = ElasticityMode.ELL_1;
                case MpccMode.elastic_ell_1_two_sided
                    obj.psi_fun_type = CFunctionType.BILINEAR_TWO_SIDED;
                    obj.relaxation_method = RelaxationMode.TWO_SIDED;
                    obj.elasticity_mode = ElasticityMode.ELL_1;
            end

            if ~obj.time_rescaling
                if obj.print_level >= 1 && obj.use_speed_of_time_variables
                    warning('nosnoc:NosnocOptions:erroneous_use_speed_of_time_variables', "use_speed_of_time_variables erroneously set to true even though we are not rescaling time, this likely means your settings are somehow faulty. Using use_speed_of_time_variables = false")
                end
                obj.use_speed_of_time_variables = 0;
            end

            if ~obj.use_speed_of_time_variables
                if obj.print_level >= 1 && obj.local_speed_of_time_variable
                    warning('nosnoc:NosnocOptions:erroneous_local_speed_of_time_variable',"local_speed_of_time_variable erroneously set to true even though we are not using speed of time variales, this likely means your settings are somehow faulty. Using local_speed_of_time_variable = false")
                end
                obj.local_speed_of_time_variable = 0;
            end

            if obj.nonlinear_sigma_rho_constraint
                obj.convex_sigma_rho_constraint = 1;
            end

            % psi function
            a = define_casadi_symbolic(obj.casadi_symbolic_mode,'a',1);
            b = define_casadi_symbolic(obj.casadi_symbolic_mode,'b',1);
            sigma = define_casadi_symbolic(obj.casadi_symbolic_mode,'sigma',1);

            switch obj.psi_fun_type
              case CFunctionType.BILINEAR
                psi_mpcc = a.*b-sigma;
                
              case CFunctionType.BILINEAR_TWO_SIDED
                psi_mpcc = [a*b-sigma,a*b+sigma];
                
              case CFunctionType.FISCHER_BURMEISTER
                psi_mpcc = a+b-sqrt(a^2+b^2+sigma^2);

              case CFunctionType.NATURAL_RESIDUAL
                psi_mpcc = 0.5*(a+b-sqrt((a-b)^2+sigma^2));

              case CFunctionType.CHEN_CHEN_KANZOW
                alpha = 0.5;
                psi_mpcc = alpha*(a+b-sqrt(a^2+b^2+sigma^2))+(1-alpha)*(a*b-sigma);

              case CFunctionType.STEFFENSON_ULBRICH
                x = a-b;
                z = x/sigma;
                y_sin = sigma*((2/pi)*sin(z*pi/2+3*pi/2)+1);
                psi_mpcc = a+b-if_else(abs(x)>=sigma,abs(x),y_sin);

              case  CFunctionType.STEFFENSON_ULBRICH_POLY
                x = a-b;
                z = x/sigma;
                y_pol = sigma*(1/8*(-z^4+6*z^2+3));
                psi_mpcc  =a+b- if_else(abs(x)>=sigma,abs(x),y_pol);

              case CFunctionType.KANZOW_SCHWARTZ
                a1 = a-sigma;
                b1 = b-sigma;
                psi_mpcc  = if_else((a1+b1)>=0,a1*b1,-0.5*(a1^2+b1^2));

              case CFunctionType.LIN_FUKUSHIMAgs
                psi_mpcc1 = [a*b-sigma];
                psi_mpcc2 = [-((a-sigma)*(b-sigma)-sigma^2)];
                psi_mpcc = vertcat(psi_mpcc1, psi_mpcc2)

              case CFunctionType.KADRANI
                psi_mpcc = (a-sigma)*(b-sigma);
            end

            obj.psi_fun = Function('psi_fun',{a,b,sigma},{psi_mpcc});

            % create Butcher tableau
            % TODO this should live somewhere else. (i.e. butcher tableu should not be in settings)
            switch obj.irk_representation
              case IrkRepresentation.integral
                [B, C, D, tau_root] = generate_butcher_tableu_integral(obj.n_s, obj.irk_scheme);
                if tau_root(end) == 1
                    right_boundary_point_explicit  = 1;
                else
                    right_boundary_point_explicit  = 0;
                end
                obj.B_irk = B;
                obj.C_irk = C;
                obj.D_irk = D;
                % also get time steps
                [~, ~, c_irk] = generate_butcher_tableu(obj.n_s,obj.irk_scheme);
                obj.c_irk = c_irk;

              case {IrkRepresentation.differential, IrkRepresentation.differential_lift_x}
                [A_irk,b_irk,c_irk,order_irk] = generate_butcher_tableu(obj.n_s,obj.irk_scheme);
                if c_irk(end) <= 1+1e-9 && c_irk(end) >= 1-1e-9
                    right_boundary_point_explicit  = 1;
                else
                    right_boundary_point_explicit  = 0;
                end
                obj.A_irk = A_irk;
                obj.b_irk = b_irk;
                obj.c_irk = c_irk;
            end
            obj.right_boundary_point_explicit = right_boundary_point_explicit;
        end

        function time_rescaling = get.time_rescaling(obj)
            time_rescaling = (obj.time_freezing && obj.impose_terminal_phyisical_time) || obj.time_optimal_problem;
        end
    end
end
