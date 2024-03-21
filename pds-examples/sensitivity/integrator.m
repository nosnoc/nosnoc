function [prob, data, opts, h] = integrator(use_fesd,N_sim, x0)
    import casadi.*
    T = 1.0;
    t_step = T/N_sim;
    %% Define (uncontrolled for now) projected system
    x = SX.sym('x', 2);
    data.x = x;
    data.lbx = [-inf;-inf];
    data.ubx = [inf;inf];
    data.x0 = [0;0.5];
    data.u = [];
    data.lbu = [];
    data.ubu = [];
    data.u0 = [];
    data.c = [x(1)+x(2)];
    data.f_x = [-4; 1];
    data.f_q = 0;%x(1)^2;
    data.f_q_T = 0;

    data.T = t_step;
    data.N_stages = 1;
    data.N_fe = 2;
    data.n_s = 2;
    data.irk_scheme = 'radau';

    h = t_step/data.N_fe;

    opts.step_eq = 'heuristic_mean';
    opts.use_fesd = use_fesd;
    opts.elastic_ell_inf = 0;

    prob = InclusionProblem(data, opts);

    prob.generate_constraints();
    f = prob.f;
    %prob.f = 0;
    %prob.g.objective(1) = {f, -inf, inf};

    default_tol = 1e-10;

    opts_casadi_nlp.ipopt.print_level = 2;
    opts_casadi_nlp.print_time = 0;
    opts_casadi_nlp.ipopt.sb = 'yes';
    opts_casadi_nlp.verbose = false;
    opts_casadi_nlp.ipopt.max_iter = 10000;
    opts_casadi_nlp.ipopt.bound_relax_factor = 0;
    %opts_casadi_nlp.ipopt.bound_relax_factor = 1e-8;
    %opts_casadi_nlp.ipopt.honor_original_bounds = 'yes';
    opts_casadi_nlp.ipopt.tol = default_tol;
    opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
    opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
    opts_casadi_nlp.ipopt.compl_inf_tol = default_tol;
    opts_casadi_nlp.ipopt.acceptable_tol = 1e-6;
    opts_casadi_nlp.ipopt.mu_strategy = 'adaptive';
    opts_casadi_nlp.ipopt.mu_oracle = 'quality-function';
    opts_casadi_nlp.ipopt.warm_start_init_point = 'yes';
    opts_casadi_nlp.ipopt.warm_start_entire_iterate = 'yes';
    opts_casadi_nlp.ipopt.linear_solver = 'ma27';
    prob.create_solver(opts_casadi_nlp);
end
