function [x_res,t_grid,model,problem_options, solver_options] = test_simple_switch(rk_representation, rk_scheme, dcs_mode, cross_comp_mode)
    import casadi.*
    problem_options = nosnoc.Options();
    integrator_options = nosnoc.integrator.Options();
    solver_options = integrator_options.fesd_solver_opts;
    model = nosnoc.model.Pss();

    problem_options.n_s = 2;
    solver_options.homotopy_update_slope = 0.1;
    problem_options.rk_scheme = rk_scheme;
    
    problem_options.cross_comp_mode = cross_comp_mode;
    problem_options.rk_representation= rk_representation;
    problem_options.dcs_mode = dcs_mode;
    solver_options.print_level = 3;
    % discretization parameters
    N_sim = 16;
    T_sim = 1.5;


    problem_options.N_sim = N_sim;
    problem_options.N_finite_elements = 2;
    problem_options.T_sim = T_sim;

    model.x0 = -1;
    x = SX.sym('x',1);
    model.x = x;
    model.c = x;
    model.S = [-1; 1];
    f_1 = [2]; f_2 = [0.2];
    model.F = [f_1 f_2];
    integrator_options.use_previous_solution_as_initial_guess = 1;
    integrator = nosnoc.Integrator(model, problem_options, integrator_options);
    [t_grid, x_res] = integrator.simulate();
end
