function x_res = integrate_system(N_stages, T_sim, u_vals)
    import casadi.*
    integrator_options = nosnoc.integrator.Options();
    solver_options = integrator_options.fesd_solver_opts;
    solver_options.complementarity_tol = 1e-10;
    solver_options.print_level = 3;
    solver_options.homotopy_steering_strategy = 'ELL_INF';
    solver_options.decreasing_s_elastic_upper_bound = true;

    [model, problem_options] = pds_ocp_dynamics(true, 1, 4, T_sim);
    problem_options.T_sim = T_sim;
    problem_options.N_sim = N_stages;
    problem_options.N_finite_elements = 4;
    
    integrator = nosnoc.Integrator(model, problem_options, solver_options);

    [t_grid,x_res,~,~] = integrator.simulate(u=u_vals, x0=model.x0)
end
