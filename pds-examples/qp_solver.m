function [result, prob] = qp_solver(Q,c,A,b,x0)
% This example implements a projected gradient method via FESD discretization
% Yes this is a little ridiculous to solve a qp by generating an MPCC and solving that with IPOPT
% which itself uses a QP solver underneath.
% Fundementally though this MPCC is (TODO(@anton) verify) linear with linear complementarity constraints which I believe should
% could be solved with a dedicated solver. Now is there any benefit here over an event capturing method? I don't think so,
% except possibly that we take "full steps" though this may actually be bad? idk :)
%
% Now also implements a rudimentary acceleration algorithm. in terms of nearly-LCP solves this actually doesn't do too badly
% if you compare to the number of iterations taken in qpoases for example, however these nearly LCP problems would need a
% custom solver to maybe make it competitive in terms of real time.
    import casadi.*
    n_x = size(x0,1);
    % TODO check sizes

    % Generate PDS corresponding to the QP
    x = SX.sym('x', n_x);
    data.x = x;
    data.lbx = -inf*ones(n_x,1);
    data.ubx = inf*ones(n_x,1);
    data.x0 = x0;
    data.u = [];
    data.lbu = [];
    data.ubu = [];
    data.u0 = [];
    data.c = A*x - b;
    data.f_x = -(Q*x + c);
    data.f_q = 0;
    data.f_q_T = 0;

    data.T = 1.0; % TODO(@anton) adjust this depending on when the homotopy fails to converge :)
    data.N_stages = 1;
    data.N_fe = 2;
    data.n_s = 1;
    data.irk_scheme = 'radau';

    % Generate problem
    prob = InclusionProblem(data, struct);

    prob.generate_constraints();

    % options
    default_tol = 1e-12;

    opts_casadi_nlp.ipopt.print_level = 2;
    opts_casadi_nlp.print_time = 0;
    opts_casadi_nlp.ipopt.sb = 'yes';
    opts_casadi_nlp.verbose = false;
    opts_casadi_nlp.ipopt.max_iter = 5000;
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

    tol = 1e-7;
    
    x_res = data.x0;
    h_res = [];
    x_curr = data.x0;
    x_prev = data.x0;
    iters = 0;
    T_max = data.T*10;
    T_curr = data.T * 1e-6;
    backoff = 0.5;
    growth = 2.0;
    success = 1;
    while true
        prob.w.x(0,0,data.n_s).init = x_curr;
        prob.w.x(0,0,data.n_s).lb = x_curr;
        prob.w.x(0,0,data.n_s).ub = x_curr;
        prob.p.T(1).init = T_curr;
        for ii=1:data.N_fe
            prob.w.h(1,ii).ub = 2*(T_curr/data.N_fe);
        end
        success = homotopy(prob);
        iters = iters+1;
        disp(['iters: ' num2str(iters) ' step size: ' num2str(T_curr)]);
        if ~success
            T_curr = T_curr*backoff;
            % TODO(@anton) can check here because incidentaly the terminal point found is sometimes actually feasible :)
            continue
        end
        T_curr = min(T_max, T_curr*growth);
        x_curr = prob.w.x(1,data.N_fe,data.n_s).res;
        lambda_curr = prob.w.lambda(1,data.N_fe,data.n_s).res;
        x_sim = prob.w.x(1,:,data.n_s).res;
        h_sim = prob.w.h(1,:).res;
        x_res = [x_res,x_sim];
        h_res = [h_res,h_sim];
        if max(abs(x_prev-x_curr)) <= tol
            result.x = x_curr;
            result.lambda = lambda_curr;
            result.traj = x_res;
            result.h_res = h_res;
            result.iters = iters;
            break
        end
        if iters > 500
            result.x = x_curr;
            result.traj = x_res;
            warning("max iters reached")
            break
        end
        x_prev = x_curr;
        prob.w.init = prob.w.res;
    end
end
