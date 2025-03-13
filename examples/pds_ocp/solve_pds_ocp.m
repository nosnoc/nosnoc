function [x_res, u_res, t_res, t_control, lambda_res, c_res] = solve_pds_ocp(use_fesd, N_stages, n_s, T, x_target)
    
    solver_options = nosnoc.solver.Options();
    solver_options.complementarity_tol = 1e-10;
    solver_options.print_level = 3;
    
    [model, problem_options] = pds_ocp_dynamics(use_fesd, N_stages, n_s, T, x_target);

    % Solve
    active_set_guess = nosnoc.activeset.Pds({[],[1]},'times', [problem_options.T*0.5,problem_options.T]);
    %active_set_guess = nosnoc.activeset.Pds({[]},'times', [problem_options.T]);
    solver_options.mpecopt.initialization_strategy = 'TakeProvidedActiveSet';
    ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
    ocp_solver.set_initial_active_set(active_set_guess);
    ocp_solver.solve('mpecopt');

    x_res = ocp_solver.get('x');
    u_res = ocp_solver.get('u');
    lambda_res = ocp_solver.get('lambda');
    c_res = ocp_solver.get('c_lift');
    t_res = ocp_solver.get_time_grid();
    t_control = ocp_solver.get_control_grid();
end
