function [x_res, u_res, t_res, t_control, lambda_res, c_res] = solve_pds_ocp(use_fesd, N_stages, n_s, T, x_target)
    
    solver_options = nosnoc.reg_homotopy.Options();
    solver_options.complementarity_tol = 1e-10;
    solver_options.print_level = 3;
    
    [model, problem_options] = pds_ocp_dynamics(use_fesd, N_stages, n_s, T, x_target);

    ccopt_rolloff = nosnoc.ccopt.Options();
    ccopt_rolloff.solver_name = 'CCOpt Rolloff';
    ccopt_rolloff.opts_madnlp.linear_solver = 'Ma27Solver';
    ccopt_rolloff.opts_madnlp.barrier.TYPE = 'MonotoneUpdate';
    ccopt_rolloff.opts_madnlp.barrier.mu_min = 1e-9;
    ccopt_rolloff.opts_ccopt.relaxation_update.TYPE = 'RolloffRelaxationUpdate';
    ccopt_rolloff.opts_ccopt.relaxation_update.rolloff_slope = 2.0;
    ccopt_rolloff.opts_ccopt.relaxation_update.rolloff_point = 1e-6;
    ccopt_rolloff.opts_ccopt.relaxation_update.sigma_min = 1e-8;
    ccopt_rolloff.opts_madnlp.tol=1e-8;
    ccopt_rolloff.opts_madnlp.print_level=3;
    ccopt_rolloff.opts_madnlp.max_iter=5000;
    ccopt_rolloff.opts_madnlp.disable_garbage_collector = true;
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
