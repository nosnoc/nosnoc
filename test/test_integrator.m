function [model,problem_options] = test_integrator(use_fesd, rk_representation, rk_scheme, dcs_mode)
    import casadi.*
    % discretization settings
    N_finite_elements = 2;
    T_sim = pi/2;
    N_sim  = 29;
    R_osc  = 1;


    fprintf('use_fesd\trk_representation\trk_scheme\tdcs_mode\n')
    fprintf('%d\t\t\t%s\t\t\t%s\t\t\t%s\n',use_fesd, rk_representation, rk_scheme, dcs_mode);
    problem_options = nosnoc.Options();
    integrator_options = nosnoc.integrator.Options();
    solver_options = integrator_options.fesd_solver_opts;
    problem_options.use_fesd = use_fesd;
    problem_options.rk_representation = rk_representation;
    problem_options.rk_scheme = rk_scheme;
    solver_options.print_level = 2;
    problem_options.n_s = 4;
    problem_options.dcs_mode = dcs_mode;
    % 'Stewart'; % 'Step;
    solver_options.complementarity_tol = 1e-9;
    problem_options.cross_comp_mode  = 3;
    solver_options.homotopy_update_rule = 'superlinear';
    solver_options.N_homotopy = 7;
    % Model
    x_star = [exp(1);0];
    x_star = [exp(T_sim-1)*cos(2*pi*(T_sim-1));-exp((T_sim-1))*sin(2*pi*(T_sim-1))];

    model = nosnoc.model.Pss();
    problem_options.N_finite_elements = N_finite_elements;
    problem_options.T_sim = T_sim;
    problem_options.N_sim = N_sim;
    omega = -2*pi;
    A1 = [1 omega;...
        -omega 1];
    A2 = [1 -omega;...
        omega 1];
    % Inital Value
    model.x0 = [exp(-1);0];
    % Variable defintion
    x1 = SX.sym('x1');
    x2 = SX.sym('x2');
    x = [x1;x2];
    c = x1^2+x2^2-1;
    model.x = x;
    model.c = c;
    model.S = [-1;1];
    f_11 = A1*x;
    f_12 = A2*x;
    F = [f_11 f_12];
    model.F = F;
    % Call integrator
    integrator = nosnoc.Integrator(model, problem_options, integrator_options);
    [t_grid, x_res] = integrator.simulate();
    % numerical error
    x_fesd = x_res(:,end);
    error_x = norm(x_fesd-x_star,"inf");
    fprintf(['Numerical error with h = %2.3f and ' char(problem_options.rk_scheme) ' with n_s = %d stages is: %5.2e: \n'],problem_options.h_sim,problem_options.n_s,error_x);
end
