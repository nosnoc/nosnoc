function [x_res,t_grid,model,problem_options, solver_options] = test_simple_cls(rk_representation, cross_comp_mode)
    import casadi.*
    %% init nosnoc settings and model
    problem_options = nosnoc.Options();
    solver_options = nosnoc.solver.Options();
    model = nosnoc.model.Cls();
    %% Simulation setings
    N_FE = 2;
    T_sim = 3;
    N_sim = 20;

    problem_options.T_sim = T_sim;
    problem_options.N_sim = N_sim;
    problem_options.N_finite_elements = N_FE;
    problem_options.rk_scheme = RKSchemes.RADAU_IIA;
    problem_options.rk_representation = rk_representation;
    problem_options.n_s = 2;
    problem_options.cross_comp_mode = cross_comp_mode;
    problem_options.dcs_mode = DcsMode.CLS;
    problem_options.no_initial_impacts = 1;

    % Initialization 
    problem_options.initial_Lambda_normal = 0;
    problem_options.initial_lambda_normal = 0;
    problem_options.initial_Y_gap = 1;
    problem_options.initial_y_gap = 1;

    %problem_options.fixed_eps_cls = 1;
    %problem_options.relax_terminal_numerical_time = ConstraintRelaxationMode.ELL_1;
    %problem_options.rho_terminal_numerical_time = 1e3;
    %problem_options.relax_fesdj_impulse = ConstraintRelaxationMode.ELL_2;
    %problem_options.rho_fesdj_impulse = 1e6;
    %problem_options.gamma_h = 0.99;

    solver_options.decreasing_s_elastic_upper_bound = 1; % elasic mode with decreasing bounds for the elstaic slacks
    solver_options.print_level = 2;


    %% model defintion
    g = 9.81;
    x0 = [0.8;0];

    q = SX.sym('q',1);
    v = SX.sym('v',1);
    model.M = 1;
    model.x = [q;v];
    model.e = 1;
    model.mu = 0;
    model.x0 = x0;
    model.f_v = -g;
    model.f_c = q;


    %% Call nosnoc integrator
    integrator = nosnoc.Integrator(model, problem_options, solver_options);
    [t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();

end
