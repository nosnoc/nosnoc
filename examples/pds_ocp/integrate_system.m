function x_res = integrate_system(N_stages, T_sim, u_vals)
    import casadi.*
    solver_options = nosnoc.solver.Options();
    solver_options.complementarity_tol = 1e-10;
    solver_options.print_level = 3;
    solver_options.homotopy_steering_strategy = 'ELL_INF';
    solver_options.decreasing_s_elastic_upper_bound = true;

    [model, problem_options] = pds_ocp_dynamics(true, 1, 4, T_sim);
    problem_options.T_sim = T_sim;
    problem_options.N_sim = N_stages;
    problem_options.N_finite_elements = 4;
    
    integrator = nosnoc.integrator.FESD(model, problem_options, solver_options);

    integrator.set_x0(model.x0);
    for ii=1:N_stages
        integrator.set('u', 'lb', {1}, u_vals(:, ii));
        integrator.set('u', 'ub', {1}, u_vals(:, ii));
        integrator.set('u', 'init', {1}, u_vals(:, ii));
        integrator.solve();
        x_next = integrator.discrete_time_problem.w.x(1,problem_options.N_finite_elements,problem_options.n_s).res;
        integrator.set_x0(x_next);
    end

    x_res = integrator.get('x');
end
