classdef NosnocOptions < handle
% TODO clean up much of the work here.
    properties
        % General
        solver_name {mustBeTextScalar} = 'nosnoc_solver'
        use_fesd(1,1) logical = 1
        casadi_symbolic_mode {mustBeMember(casadi_symbolic_mode,{'SX', 'MX'})} = 'SX' % TODO enum

        % Interface settings
        real_time_plot(1,1) logical = 0

        % IRK and FESD Settings
        n_s(1,1) {mustBeInteger, mustBePositive} = 2
        irk_scheme = 'radau' % TODO enum
        lift_irk_differential(1,1) logical = 1
        cross_comp_mode {mustBeInRange(cross_comp_mode, 1, 12)}= 3 
        gamma_h = 1
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
        n_depth_step_lifting(1,1) {mustBeInteger} = 2

        % TODO: make this static
        list_of_all_rk_schemes = {'radau','legendre','Radau-IIA','Gauss-Legendre','Radau-I','Radau-IA',...
                                  'Lobatto-III','Lobatto-IIIA','Lobatto-IIIB','Lobatto-IIIC',...
                                  'Explicit-RK'}
        %General NLP/OCP Settings
        g_path_constraint(1,1) logical = 0 % is nonlinear path constraint present (by default evaluated only on control grid points)
        g_comp_path_constraint(1,1) logical = 0
        g_path_at_fe(1,1) logical = 0 % evaluate nonlinear path constraint at every finte element boundary
        g_path_at_stg(1,1) logical = 0 % evaluate nonlinear path constraint at every stage 

        x_box_at_fe(1,1) logical = 1 % evaluate box constraint for diff states at every finite element boundary point
        x_box_at_stg(1,1) logical = 1 % evaluate box constraint for diff states at every stage point. (is set to zero per default in differential irk mode, as it becomes a linear instead of box constraint)

        terminal_constraint(1,1) logical = 0
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
        N_homotopy(1,1) {mustBeInteger, mustBePositive} = ceil(abs(log(1e-9)/log(0.1)))
        s_elastic_max(1,1) double {mustBeReal, mustBePositive} = 1e1
        s_elastic_min(1,1) double {mustBeReal, mustBeNonnegative} = 0
        s_elastic_0(1,1) double {mustBeReal, mustBePositive} = 1

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
        h_fixed_iterations(1,1) logical = 0
        h_fixed_max_iter(1,1) logical = 1 % number of iterations that are done with fixed h in the homotopy loop
        h_fixed_change_sigma(1,1) logical = 1 % if this is on, do not update sigma and just solve on nlp with fixed h.
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
        time_freezing(1,1) logical  = 0
        time_freezing_inelastic(1,1) logical = 0
        time_rescaling(1,1) logical = 0
        % for time optimal problems and equidistant control grids in physical time
        use_speed_of_time_variables(1,1) logical = 1
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

        friction_exists(1,1) logical = 0
        % exerimentla expert otpions
        nonsmooth_switching_fun(1,1) logical = 0 % experimental: use c = max(c1,c2) insetad of c = [c1c2]  

        % expert mode: stabilize auxiliary dynamics in \nabla f_c(q) direction
        stabilizing_q_dynamics(1,1) logical = 0
        kappa_stabilizing_q_dynamics(1,1) double {mustBeReal, mustBePositive} = 1e-5
        % Verbose
        print_level = 3

        % IPOPT Settings
        tol_ipopt(1,1) double {mustBeReal, mustBePositive} = 1e-16
        opts_ipopt

        % Relaxation of terminal constraint
        relax_terminal_constraint(1,1) logical = 0 %  0  - hard constraint, 1 - ell_1 , 2  - ell_2 , 3 - ell_inf
        relax_terminal_constraint_from_above(1,1) logical = 0 
        rho_terminal(1,1) double {mustBeReal, mustBePositive} = 1e2
        relax_terminal_constraint_homotopy(1,1) logical = 0 % terminal penalty is governed by homotopy parameter

        % Integrator Specific 
        use_previous_solution_as_initial_guess(1,1) logical = 0
        simulation_problem(1,1) logical = 0

        % Misc
        there_exist_free_x0(1,1) logical = 0
        time_freezing_model_exists(1,1) logical = 0


        % All NLP parameters
        T_val(1,1) double {mustBeReal, mustBePositive} = 1
        p_val

        % Butcher Tableu
        A_irk double
        B_irk double
        b_irk double
        C_irk double
        D_irk double

        right_boundary_point_explicit(1,1) logical % TODO this shoud live in model probably
    end

    methods
        function obj = NosnocOptions()
            
            obj.opts_ipopt.ipopt.print_level = 0;
            obj.opts_ipopt.print_time = 0;
            obj.opts_ipopt.ipopt.sb = 'yes';
            obj.opts_ipopt.verbose = false;
            obj.opts_ipopt.ipopt.max_iter = 500;
            % opts_ipopt.ipopt.print_level = 5
            obj.opts_ipopt.ipopt.bound_relax_factor = 0;
            obj.opts_ipopt.ipopt.tol = 1e-16;
            obj.opts_ipopt.ipopt.dual_inf_tol = 1e-16;
            obj.opts_ipopt.ipopt.dual_inf_tol = 1e-16;
            obj.opts_ipopt.ipopt.compl_inf_tol = 1e-16;
            obj.opts_ipopt.ipopt.mu_strategy = 'adaptive';
            obj.opts_ipopt.ipopt.mu_oracle = 'quality-function';

            obj.p_val = [obj.sigma_0,obj.rho_sot,obj.rho_h,obj.rho_terminal,obj.T_val]
        end

        function [] = create_butcher_tableu(obj, model)
            switch obj.irk_representation
              case 'integral'
                [B, C, D, tau_root] = generate_butcher_tableu_integral(model.dimensions.n_s, obj.irk_scheme);
                if tau_root(end) == 1
                    right_boundary_point_explicit  = 1;
                else
                    right_boundary_point_explicit  = 0;
                end
                obj.B_irk = B;
                obj.C_irk = C;
                obj.D_irk = D;
              case {'differential', 'differential_lift_x'}
                [A_irk,b_irk,c_irk,order_irk] = generate_butcher_tableu(model.dimensions.n_s,obj.irk_scheme);
                if c_irk(end) <= 1+1e-9 && c_irk(end) >= 1-1e-9
                    right_boundary_point_explicit  = 1;
                else
                    right_boundary_point_explicit  = 0;
                end
                obj.A_irk = A_irk;
                obj.b_irk = b_irk;
              otherwise
                error('Choose irk_representation either: ''integral'' or ''differential''')
            end
            obj.right_boundary_point_explicit = right_boundary_point_explicit;
        end
        
    end
    
end
