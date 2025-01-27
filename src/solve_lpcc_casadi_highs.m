function  [f_lin_opt, d_lpcc, y_lpcc, lpcc_solution_exists, lpcc_cpu_time] = solve_lpcc_casadi_highs(x_k,...
            dims,...
            d_lpcc, y_lpcc, Delta_TR_k_l,...
            lbx, ubx,...
            nabla_f_k, A_k, a_k, B_k, b_k, B_comp_k, b_comp_k)
    import casadi.*
    A_highs = [A_k,zeros(dims.n_eq,dims.n_comp);
        B_k,zeros(dims.n_ineq,dims.n_comp);
        B_comp_k];
    
    b_highs = [-a_k; -b_k; b_comp_k];

    c_highs = [nabla_f_k;zeros(dims.n_comp,1)];

    discrete = [false(dims.n_primal, 1); true(dims.n_comp, 1)];

    lb = lbx;
    lb(dims.ind_x0) = max([lbx(dims.ind_x0)-x_k(dims.ind_x0),-Delta_TR_k_l*ones(size(dims.ind_x0))]');
    lb(dims.ind_x1) = max([lbx(dims.ind_x1)-x_k(dims.ind_x1),-Delta_TR_k_l*ones(dims.n_comp,1),-x_k(dims.ind_x1)],[],2);
    lb(dims.ind_x2) = max([lbx(dims.ind_x2)-x_k(dims.ind_x2),-Delta_TR_k_l*ones(dims.n_comp,1),-x_k(dims.ind_x2)],[],2);
    lb = [lb;-inf*ones(dims.n_comp,1)];
    % ub
    ub = [min(ubx-x_k,Delta_TR_k_l);inf*ones(dims.n_comp,1)];

    lba = b_highs;
    uba = [-a_k; inf(dims.n_ineq + 2*dims.n_comp,1)];
    lp.a = DM(A_highs).sparsity();
    lpopts = struct;
    lpopts.discrete = discrete;
    %lpopts.highs.simplex_strategy = 4;
    lpopts.highs.log_to_console = false;
    lpopts.error_on_fail = false;
    lpsol = conic('lp', 'highs', lp, lpopts);

    %% solve
    r = lpsol('g', c_highs, 'a', A_highs, 'lbx', lb, 'ubx', ub, 'lba', lba, 'uba', uba);

    if isequal(lpsol.stats.return_status,'Optimal')
        lpcc_solution_exists = 1;
        d_lpcc = full(r.x(1:dims.n_primal));
        y_lpcc = full(r.x(end-dims.n_comp+1:end));
        f_lin_opt = full(r.cost);
    else
        lpcc_solution_exists = 0;
        f_lin_opt = nan;
    end
    lpcc_cpu_time = lpsol.stats.t_proc_solver;
end
