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
classdef Options < handle
% The main options class for `nosnoc`. It contains the options for simulation and optimal control reformulations but *does not*
% contain settings for the MPCC solver.
    properties
        % boolean: If true the FESD discretization is used, otherwise a direct time-stepping discretization is used.
        use_fesd(1,1) logical = 1;

        % string: Which casadi symbolics to use. Can either be `'casadi.SX'` or `'casadi.MX'.`
        casadi_symbolic_mode {mustBeMember(casadi_symbolic_mode,{'casadi.SX', 'casadi.MX'})} = 'casadi.SX';

        T_sim
        N_sim
        h_sim

        h % double: Control stage Step size.
        h_k % double: Finite element step size.
        T % double: Terminal time.

        N_stages(1,1) {mustBeInteger, mustBePositive} = 1; % int: Number of control stages.

        % int: Number of finite elements in each control stage. This can either be a scalar value
        % in which case it is transformed into a vector for that value when :meth:`preprocess` is called.
        % Alternatively you can pass a vector of size :attr:`N_stages`.
        N_finite_elements {mustBeInteger, mustBePositive} = 2;

        n_s(1,1) {mustBeInteger} = 2 % int: Number of Stages in the Runge-Kutta scheme.

        % RKSchemes: Which Runge-Kutta scheme family to use.
        %
        % See Also:
        %    `RKSchemes` for more details as to the how to choose a Runge-Kutta Scheme and
        %    for differences between them.
        rk_scheme(1,1) RKSchemes = RKSchemes.RADAU_IIA

        % RKRepresentation: Which representation of Runge-Kutta discretization to use.
        %
        % See Also:
        %     `RKRepresentation` for a description of the representations.
        rk_representation RKRepresentation = RKRepresentation.integral;

        % CrossCompMode: Which cross complementarity mode to use.
        %
        % See Also:
        %     `CrossCompMode` for a description of the representations.
        cross_comp_mode(1,1) CrossCompMode = CrossCompMode.FE_STAGE

        % double: Fraction in the range $\gamma_h \in [0,1]$ by which the step size is relaxed:
        % $$(1-\gamma_h) h_0\le h \le (1+\gamma_h) h_0$$
        gamma_h(1,1) double {mustBeReal} = 1
        dcs_mode DcsMode = DcsMode.Stewart % DcsMode: Which DCS to reformulate the problem into.

        lift_complementarities(1,1) logical = 0 % boolean: Whether complementarities are lifted. TODO(@anton) should this still live in MPCC generation?
        lower_bound_comp_lift(1,1) logical = 0 % boolean: If true we add additional lower bounds to the lifted variables.

        %--------------------- Initial Values ---------------------%
        
        initial_alpha(1,1) double {mustBeReal, mustBeFinite} = 0.5 % double: Initial value for $\alpha$ in the Heaviside step reformulation.
        initial_lambda_n(1,1) double {mustBeReal, mustBeFinite} = 0.5 % double: Initial value for $\lambda_n$ in the Heaviside step reformulation.
        initial_lambda_p(1,1) double {mustBeReal, mustBeFinite} = 0.5 % double: Initial value for $\lambda_p$ in the Heaviside step reformulation.
        initial_beta_lift(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for $\beta$ when lifting is enabled in the Heaviside step reformulation.
        initial_theta_step(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for $\theta$ when lifting is enabled in the Heaviside step reformulation.
        initial_lambda_gcs(1,1) double {mustBeReal, mustBeFinite} = 0 % double: Initial value for $\lambda$ in the Gradient Comlementarity System.
        
        initial_Lambda_normal(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for $\Lambda_n$ in FESD-J reformulation.
        initial_P_vn(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for positive normal velocity slack in FESD-J reformulation impulse calculation.
        initial_N_vn(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for negative normal velocity slack in FESD-J reformulation impulse calculation.
        initial_Y_gap(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for gap function in FESD-J reformulation impulse calculation.
        initial_Lambda_tangent(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for $\Lambda_t$ in FESD-J reformulation impulse calculation.
        initial_Gamma_d(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for $\Gamma_d$ in FESD-J reformulation impulse calculation.
        initial_Beta_d(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for $\Beta_d$ in FESD-J reformulation impulse calculation.
        initial_Delta_d(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for $\Delta_d$ in FESD-J reformulation impulse calculation.
        initial_Gamma(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for $\Gamma$ in FESD-J reformulation impulse calculation.
        initial_Beta(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for $\Beta$ in FESD-J reformulation impulse calculation.
        initial_P_vt(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for positive tangential velocity slack in FESD-J reformulation impulse calculation.
        initial_N_vt(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for negative tangential velocity slack in FESD-J reformulation impulse calculation.
        initial_Alpha_vt(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value fo tangential velocity step function in FESD-J reformulation impulse calculation.

        initial_lambda_normal(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for $\lambda_n$ in FESD-J reformulation.
        initial_p_vn(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for positive normal velocity slack in FESD-J reformulation.
        initial_n_vn(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for negative normal velocity slack in FESD-J reformulation.
        initial_y_gap(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for gap function in FESD-J reformulation.
        initial_lambda_tangent(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for $\lambda_t$ in FESD-J reformulation.
        initial_gamma_d(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for $\gamma_d$ in FESD-J reformulation.
        initial_beta_d(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for $\beta_d$ in FESD-J reformulation.
        initial_delta_d(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for $\delta_d$ in FESD-J reformulation.
        initial_gamma(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for $\gamma$ in FESD-J reformulation.
        initial_beta(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for $\beta$ in FESD-J reformulation.
        initial_p_vt(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for positive tangential velocity slack in FESD-J reformulation.
        initial_n_vt(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value for negative tangential velocity slack in FESD-J reformulation.
        initial_alpha_vt(1,1) double {mustBeReal, mustBeFinite} = 1 % double: Initial value fo tangential velocity step function in FESD-J reformulation.
        
        %--------------------- End Initial Values ---------------------%

        %--------------------- Max Values ---------------------%
        
        ub_lambda_gcs(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for $\lambda$ in the Gradient Comlementarity System.
        
        ub_Lambda_normal(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for $\Lambda_n$ in FESD-J reformulation.
        ub_P_vn(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for positive normal velocity slack in FESD-J reformulation impulse calculation.
        ub_N_vn(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for negative normal velocity slack in FESD-J reformulation impulse calculation.
        ub_Y_gap(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for gap function in FESD-J reformulation impulse calculation.
        ub_Lambda_tangent(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for $\Lambda_t$ in FESD-J reformulation impulse calculation.
        ub_Gamma_d(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for $\Gamma_d$ in FESD-J reformulation impulse calculation.
        ub_Beta_d(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for $\Beta_d$ in FESD-J reformulation impulse calculation.
        ub_Delta_d(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for $\Delta_d$ in FESD-J reformulation impulse calculation.
        ub_Gamma(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for $\Gamma$ in FESD-J reformulation impulse calculation.
        ub_Beta(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for $\Beta$ in FESD-J reformulation impulse calculation.
        ub_P_vt(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for positive tangential velocity slack in FESD-J reformulation impulse calculation.
        ub_N_vt(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for negative tangential velocity slack in FESD-J reformulation impulse calculation.
        ub_Alpha_vt(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value fo tangential velocity step function in FESD-J reformulation impulse calculation.

        ub_lambda_normal(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for $\lambda_n$ in FESD-J reformulation.
        ub_p_vn(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for positive normal velocity slack in FESD-J reformulation.
        ub_n_vn(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for negative normal velocity slack in FESD-J reformulation.
        ub_y_gap(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for gap function in FESD-J reformulation.
        ub_lambda_tangent(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for $\lambda_t$ in FESD-J reformulation.
        ub_gamma_d(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for $\gamma_d$ in FESD-J reformulation.
        ub_beta_d(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for $\beta_d$ in FESD-J reformulation.
        ub_delta_d(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for $\delta_d$ in FESD-J reformulation.
        ub_gamma(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for $\gamma$ in FESD-J reformulation.
        ub_beta(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for $\beta$ in FESD-J reformulation.
        ub_p_vt(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for positive tangential velocity slack in FESD-J reformulation.
        ub_n_vt(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value for negative tangential velocity slack in FESD-J reformulation.
        ub_alpha_vt(1,1) double {mustBeReal, mustBePositive} = inf % double: Max value fo tangential velocity step function in FESD-J reformulation.

        %--------------------- End Max Values ---------------------%
       
        % boolean: If true then the convex multiplier expressions are lifted in the Heaviside step reformulation.
        %
        % Warning:
        %     This is not currently implemented for generic Heaviside step DCS.
        pss_lift_step_functions(1,1) logical = 0 
        n_depth_step_lifting(1,1) {mustBeInteger, mustBeGreaterThanOrEqual(n_depth_step_lifting, 2)} = 2 % int: Depth to which the Heaviside step convex multipliers are lifted.

        gcs_lift_gap_functions(1,1) logical = 1 % boolean: If true the step functions $c(x)$ are lifted in the gradient complementarity system reformulation.

        linear_complemtarity_M(1,1) double {mustBeReal, mustBePositive} = 1000 % double: $M$ multiplier used in the linear complementarity step equilibration fromulation. Larger values alleviate infeasibility for smaller step sizes.

        g_path_at_fe(1,1) logical = 0 % boolean: If true we evaluate nonlinear path constraint at every finte element boundary.
        g_path_at_stg(1,1) logical = 0 % boolean: If true evaluate nonlinear path constraint at every stage.

        x_box_at_fe(1,1) logical = 1 % boolean: If true we evaluate box constraint for diff states at every finite element boundary point.

        % boolean: If true we evaluate box constraint for diff states at every stage point.
        %
        % Note:
        %    This is set to zero per default in differential rk mode, as it becomes a linear instead of box constraint.
        x_box_at_stg(1,1) logical = 1
        time_optimal_problem(1,1) = 0 % boolean: If true for an OCP we automatically reformulate the problem to be time optimal.

        rho_h(1,1) double {mustBeNonnegative} = 1 % double: Weight used in heuristic or relaxed step equilibration modes.

        % StepEquilibrationMode: Which step equilibration mode to use.
        %
        % See Also:
        %     `StepEquilibrationMode` for more details on how each mode works.
        step_equilibration(1,1) StepEquilibrationMode = StepEquilibrationMode.heuristic_mean
        step_equilibration_sigma(1,1) double {mustBePositive} = 0.1 % double: Slope at zero for the sigmoid used to rescale the indicator function, nu_ki_rescaled = tanh(nu_ki/step_equilibration_sigma).

        equidistant_control_grid(1,1) logical = 1 % boolean: If true each control stage is fixed length.

        time_freezing(1,1) logical = 0 % boolean: Use a time freezing reformulation for the given model.
        time_freezing_inelastic(1,1) logical = 0 % boolean: Use the specailized time freezing reformulation for systems with inelastic collisions and friction. 
        
        use_speed_of_time_variables(1,1) logical = 0 % boolean: If true speed of time variables are used for the time freezing reformulation or time optimal problem
        local_speed_of_time_variable(1,1) logical = 0 % boolean: If true then each control stage has a speed of time variable. Otherwise a single speed of time variable is used.
        stagewise_clock_constraint(1,1) logical = 1 % boolean: If true the control grid is fixed with constraints for each control stage.
        impose_terminal_phyisical_time(1,1) logical = 1 % boolean: If true the terminal physical time in a time freezing system is constrained to be exactly the desired horizon length.
        s_sot0(1,1) double {mustBePositive} = 1 % double: Initial value for speed of time variables.
        s_sot_max(1,1) double {mustBePositive} = 25 % double: Maximum for speed of time variables.
        s_sot_min(1,1) double {mustBePositive} = 1 % double: Minimum for speed of time variables.
        S_sot_nominal(1,1) double {mustBePositive} = 1 % double: Nominal speed of time used for regularizing the speed of time variables.
        rho_sot(1,1) double {mustBeReal, mustBeNonnegative} = 0 % double: Weight used for the speed of time regularization.

        T_final_max(1,1) double {mustBePositive} = 1e2 % double: Maximum final time for a time optimal problem.
        T_final_min(1,1) double {mustBeReal, mustBeNonnegative} = 0 % double: Minimum final time for a time optimal problem.


        time_freezing_reduced_model(1,1) logical = 0 % boolean: Analytic reduction of lifter formulation, less algebraic variables (experimental). TODO(@armin) What was this supposed to be?
        time_freezing_hysteresis(1,1) logical = 0 
        time_freezing_nonlinear_friction_cone(1,1) logical = 1 % boolean: If true we use the nonlinear friction cone, otherwise use polyhedral l_inf approximation.

        time_freezing_quadrature_state(1,1) logical = 0 % boolean: If true make a nonsmooth quadrature state to integrate only if physical time is running.
        time_freezing_lift_forces(1,1) logical = 0 % If true replace $\dot{v} = M(q)^{-1}f(q,v,u)$ by $dot{v} = z,  M(q)z - f(q,v,u) = 0$.

        % boolean: Experimental, use $c = \max(c1,c2)$ insetad of $c = c_1c_2$.
        % This is used to reduce the number of switching functions needed to generate the T shaped intersections
        % in inelastic time freezing reformulation.
        time_freezing_nonsmooth_switching_fun(1,1) logical = 0

        % boolean: Stabilize auxiliary dynamics in \nabla f_c(q) direction in the style of Baumgartner stabilization.
        stabilizing_q_dynamics(1,1) logical = 0

        % double: Constant used for stabilizing auxiliary dynamics in \nabla f_c(q) direction.
        kappa_stabilizing_q_dynamics(1,1) double {mustBePositive} = 1e-5

        % int: Level of verbosity that the `nosnoc` reformulator uses.
        %
        % Todo:
        %    @anton, @armin document this better.
        print_level = 3

        % FrictionModel: Which Friction model to use for the Complementarity Lagrangian System.
        %
        % Default: :mat:class:`FrictionModel.Conic`
        %
        % See Also:
        %     `FrictionModel` for more details as to the differences between the friction models.
        friction_model (1,1) FrictionModel = FrictionModel.Conic;

        % ConicModelSwitchHandling: Which velocity switch handling mode to use when using the Conic friction model
        %
        % See Also:
        %     `ConicModelSwitchHandling` for more details as to the differences between the switch handling modes.
        conic_model_switch_handling (1,1) ConicModelSwitchHandling = ConicModelSwitchHandling.Abs;
        
        kappa_friction_reg (1,1) double {mustBeReal, mustBeNonnegative} = 0; % double: Regularization term in friction equations to avoid large multipliers if no contact happens.
        
        lift_velocity_state(1,1) logical = 0; % boolean: If true define auxliary algebraic vairable, $dot{v} = z_v$, to avoid symbolic inversion of the inertia matrix.
        eps_cls double = 1e-3 % double: enforce $f_c$ at Euler step with h * eps_cls
        fixed_eps_cls logical = false % boolean: use fixed step eps_cls instead of a multiple of h.

        % double: The constant radius of relaxation for the friction force which enforces a nonempty interior around zero velocity
        %
        % See Also:
        %     More details can be found in :cite:p:`Nurkanovic2023a`
        eps_t double = 1e-7
        
        % ConstraintRelaxationMode: What (if any) relaxation to apply to the terminal constraints.
        %
        % See Also:
        %    `ConstraintRelaxationMode` for a detailed description of the available relaxation modes.
        relax_terminal_constraint(1,1) ConstraintRelaxationMode = ConstraintRelaxationMode.NONE; 
        relax_terminal_constraint_from_above(1,1) logical = 0; % boolean: If true we only relax the upper bound of the terminal constraint. TODO(@armin) do we still want this. I have never seen it be useful.
        rho_terminal(1,1) double {mustBePositive} = 1e2; % double: Weight used to penalize terminal constraint violation.

        % boolean: If True the terminal constraint violation penalty is governed by homotopy parameter.
        %
        % Warning:
        %     This option is currently unimplemented.
        relax_terminal_constraint_homotopy(1,1) logical = 0; 

        % ConstraintRelaxationMode: What (if any) relaxation to apply to the terminal/or stage numerical time constraints.
        %
        % See Also:
        %    `ConstraintRelaxationMode` for a detailed description of the available relaxation modes.
        relax_terminal_numerical_time(1,1) ConstraintRelaxationMode = ConstraintRelaxationMode.NONE;
        rho_terminal_numerical_time(1,1) double {mustBeNonnegative} = 1e2 % double: Weight used to penalize terminal numerical time violation.
        
        % boolean: If True the terminal numerical time constraint violation penalty is governed by homotopy parameter
        %
        % Warning:
        %     This option is currently unimplemented
        relax_terminal_numerical_time_homotopy (1,1) logical = 0; % us the homotopy parameter for the penalty.

        % ConstraintRelaxationMode: What (if any) relaxation to apply to the terminal/or stage phyical time constraints.
        %
        % See Also:
        %    `ConstraintRelaxationMode` for a detailed description of the available relaxation modes.
        relax_terminal_physical_time(1,1) ConstraintRelaxationMode = ConstraintRelaxationMode.NONE; % instead of imposing $t(T) = T$, add it as $\ell_1$ penalty term.
        rho_terminal_physical_time(1,1) double {mustBeNonnegative} = 1e2 % double: Weight used to penalize terminal physical time violation.

        % boolean: If True the terminal physical time constraint violation penalty is governed by homotopy parameter.
        %
        % Warning:
        %     This option is currently unimplemented.
        relax_terminal_physical_time_homotopy (1,1) logical = 0;

        % boolean: If false then the Lagrange term is integrated correctly, otherwise we only evaluate it at the
        % ends of control stages. Setting this to true allows us to do parameter estimation with a nonlinear cost function.
        % This is useful to set to true when implementing a maximum liklihood estimator as in combination with
        % an equidistant grid it allows for fixed time grid for measurements.
        euler_cost_integration(1,1) logical = 0 

        no_initial_impacts(1,1) logical = 0 % boolan: If true we disallow impulsive contacts at the beginning of the first control stage. 

        use_previous_solution_as_initial_guess(1,1) logical = 0 % boolean: When simulating use the previous step as an initial guess for the current one.

        has_clock_state(1,1) logical = 0
        
        T_val(1,1) double {mustBePositive} = 1
        p_val

        % Time Freezing constants
        a_n(1,1) double {mustBePositive} = 100;
        k_aux(1,1) double {mustBePositive} = 10;
        time_freezing_Heaviside_lifting(1,1) logical = true; % boolean: Exploit the time-freezing PSS structure for tailored lifting in Heaviside reformulation, and drastically reduce the number of  algebraic variables.
        
        % Butcher Tableu
        A_rk double
        B_rk double
        b_rk double
        C_rk double
        D_rk double
        c_rk double

        right_boundary_point_explicit(1,1) logical

        % experimental:
        %---------------------------------------------------------------------%

        use_numerical_clock_state(1,1) logical = false % logical: instead of sum of $h$ being used for equidistant control steps use a simple integrated state.
    end        


    properties(Dependent)
        time_rescaling
    end

    methods
        function obj = Options()
            check_matlab_requirement()
            obj.p_val = [obj.rho_sot, obj.rho_h, obj.rho_terminal, obj.T_val];
        end

        function [] = preprocess(obj)
            import casadi.*

            % time grid
            % TODO: merge T_sim and T?
            if ~isempty(obj.N_sim) && ~isempty(obj.T_sim)
                obj.T = obj.T_sim/obj.N_sim;
                obj.h_sim = obj.T_sim/(obj.N_sim*obj.N_stages*obj.N_finite_elements);
                if obj.print_level >= 2 && exist("h_sim")
                    fprintf('Info: N_sim is given, so the h_sim provided by the user is overwritten.\n')
                end
            elseif ~isempty(obj.N_sim) || ~isempty(obj.T_sim)
                error('Provide both N_sim and T_sim for the integration.')
            end

            if numel(obj.T) ~= 1 && ~obj.time_optimal_problem
                error('terminal numerical time T must be provided if time_optimal_problem is False.');
            elseif numel(obj.T) == 0 && ~obj.time_optimal_problem
                obj.T = 1;
            elseif numel(obj.T) ~= 1
                error('terminal time T must be a positive scalar.');
            end
            obj.h = obj.T/obj.N_stages;

            % check irk scheme compatibility
            if ismember(obj.rk_scheme, RKSchemes.differential_only)
                if obj.print_level >=1
                    fprintf(['Info: The user provided RK scheme: ' char(obj.rk_scheme) ' is only available in the differential representation.\n']);
                end
                obj.rk_representation = RKRepresentation.differential;
            end
            if obj.n_s < 1 || obj.n_s > 9
                error("n_s must be in [1, 9]");
            end


            if obj.gamma_h < 0 || obj.gamma_h > 1
                error("gamma_h must be in [0, 1]");
            end

            % check impacts mode
            if obj.dcs_mode == DcsMode.CLS && obj.time_freezing == 1
                if obj.print_level >=1
                    fprintf(['nosnoc: User uses dcs_mode.CLS and time_freezing = true at the same time. Defaulting to time freezing as it is currently more competitive\n']);
                end
                obj.time_freezing = 1;
                obj.dcs_mode = DcsMode.Heaviside;
            end

            if obj.time_freezing
                obj.use_speed_of_time_variables = 1;
                obj.local_speed_of_time_variable = 1;
            end

            % N finite elements shape \
            if isscalar(obj.N_finite_elements)
                obj.N_finite_elements = obj.N_finite_elements*ones(obj.N_stages,1);
            else
                if length(obj.N_finite_elements) == obj.N_stages
                    obj.N_finite_elements = obj.N_finite_elements;
                else
                    error('settings.N_finite_elements must be length 1 or N_stages');
                end
            end

            obj.h_k = obj.h./obj.N_finite_elements;

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

            % create Butcher tableau
            % TODO this should live somewhere else. (i.e. butcher tableu should not be in settings)
            switch obj.rk_representation
              case RKRepresentation.integral
                [B, C, D, tau_root] = generate_butcher_tableu_integral(obj.n_s, obj.rk_scheme);
                if tau_root(end) == 1
                    right_boundary_point_explicit  = 1;
                else
                    right_boundary_point_explicit  = 0;
                end
                obj.B_rk = B;
                obj.C_rk = C;
                obj.D_rk = D;
                % also get time steps
                [~, ~, c_rk] = generate_butcher_tableu(obj.n_s,obj.rk_scheme);
                obj.c_rk = c_rk;

              case {RKRepresentation.differential, RKRepresentation.differential_lift_x}
                [A_rk,b_rk,c_rk,order_irk] = generate_butcher_tableu(obj.n_s,obj.rk_scheme);
                if c_rk(end) <= 1+1e-9 && c_rk(end) >= 1-1e-9
                    right_boundary_point_explicit  = 1;
                else
                    right_boundary_point_explicit  = 0;
                end
                obj.A_rk = A_rk;
                obj.b_rk = b_rk;
                obj.c_rk = c_rk;
            end
            if obj.rk_representation == RKRepresentation.differential
                obj.x_box_at_stg = 0;
            end
            obj.right_boundary_point_explicit = right_boundary_point_explicit;

            if obj.N_stages == 1
                obj.stagewise_clock_constraint = false;
            end
        end

        function time_rescaling = get.time_rescaling(obj)
            time_rescaling = (obj.time_freezing && obj.impose_terminal_phyisical_time) || obj.time_optimal_problem;
        end
    end
end
