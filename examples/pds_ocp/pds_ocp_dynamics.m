function [model, problem_options] = pds_ocp_dynamics(use_fesd, N_stages, n_s, T, x_target)
    import casadi.*
    problem_options = nosnoc.Options();
    problem_options.rk_scheme = RKSchemes.RADAU_IIA;
    problem_options.rk_representation = RKRepresentation.integral;
    problem_options.cross_comp_mode = CrossCompMode.FE_FE;
    problem_options.N_finite_elements = 2;
    problem_options.n_s = n_s;
    problem_options.N_stages = N_stages;
    problem_options.T = T;
    problem_options.rho_h = 1e-10;
    problem_options.use_fesd = use_fesd;
    problem_options.gcs_lift_gap_functions = true;


    model = nosnoc.model.Pds();
    x = SX.sym('x',2);
    model.x = [x];
    model.lbx = [-inf;-inf];
    model.ubx = [inf;inf];
    x0 =[1; 5];
    model.x0 = [x0];
    u = SX.sym('u1', 2);;
    model.u = [u];
    model.lbu = [-1;-1];
    model.ubu = [1;1];
    model.u0 = [0;0];
    P = [1/4, 0;
        0, 1/16];
    model.c = [x'*P*x - 1];
    model.f_x_unconstrained = [-0.2*(x(1)+1)^2;-0.4*(x(2)+3)] + u;

    if exist('x_target')
        % costs
        R = diag([1e0;1e0]);
        Q_T = diag([0;0]);
        model.f_q = u'*R*u;
        model.f_q_T = 0.5*(x-x_target)'*Q_T*(x-x_target);
        model.g_terminal = x-x_target;
    end
end
