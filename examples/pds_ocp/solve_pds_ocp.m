function [x_res, u_res, h_k] = solve_pds_ocp(use_fesd, N_stages, n_s, T, x_target)
    
    solver_options = nosnoc.solver.Options();
    solver_options.complementarity_tol = 1e-10;
    solver_options.print_level = 3;
    
    [model, problem_options] = pds_ocp_dynamics(use_fesd, N_stages, n_s, T, x_target);

    % Solve
    ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
    ocp_solver.solve();

    h_k = problem_options.h_k(end);
    x_res = ocp_solver.get('x');
    u_res = ocp_solver.get('u');
end
