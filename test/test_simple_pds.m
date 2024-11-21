function [x_res,t_grid,model,problem_options, solver_options] = test_simple_pds(rk_representation, cross_comp_mode)
    import casadi.*
    N_sim = 31;
    T_sim = 11*pi/12 + sqrt(3);
    N_finite_elements = 2;
    %% nosnoc settings
    problem_options = nosnoc.Options();
    integrator_options = nosnoc.integrator.Options();
    solver_options = integrator_options.fesd_solver_opts;

    problem_options.n_s = 3;
    problem_options.rk_scheme = RKSchemes.RADAU_IIA;
    problem_options.rk_representation = rk_representation;
    problem_options.cross_comp_mode = cross_comp_mode;
    problem_options.N_sim = N_sim;
    problem_options.N_finite_elements = N_finite_elements;
    problem_options.T_sim = T_sim;
    problem_options.gcs_lift_gap_functions = true;
    %solver_options.homotopy_steering_strategy = 'ELL_INF';
    solver_options.complementarity_tol = 1e-10;
    solver_options.print_level = 2;

    x = SX.sym('x',2);
    model = nosnoc.model.Pds();
    model.x0 = [sqrt(2);sqrt(2)];
    model.x = x;
    model.c = x(2)+1;
    model.f_x_unconstrained = [x(2);-x(1)];

    integrator = nosnoc.Integrator(model, problem_options, integrator_options);
    [t_grid, x_res] = integrator.simulate();
end
