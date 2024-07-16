function [x_res, u_res, t_res, t_control, lambda_res, c_res] = solve_pds_ocp(use_fesd, N_stages, n_s, T, x_target)
    
    solver_options = nosnoc.solver.Options();
    solver_options.complementarity_tol = 1e-10;
    solver_options.print_level = 3;
    
    [model, problem_options] = pds_ocp_dynamics(use_fesd, N_stages, n_s, T, x_target);

    % Solve
    ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
    ocp_solver.solve();

    x_res = ocp_solver.get('x');
    u_res = ocp_solver.get('u');
    lambda_res = ocp_solver.get('lambda');
    c_res = ocp_solver.get('c_lift');
    t_res = ocp_solver.get_time_grid();
    t_control = ocp_solver.get_control_grid();
end
