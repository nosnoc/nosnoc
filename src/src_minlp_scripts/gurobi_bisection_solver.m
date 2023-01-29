%% Gurobi Bisection solver
miqp_casadi_functions
sucess = 0;
CPU_time_gurobi = 0;
T_opt = T_ub_k;
lb_k = [T_lb_k];
ub_k = [T_ub_k];
T_val_k = T_ub_k;
diff_ub_lb = ub_k(end)-lb_k(end);
iter = 0;
%     for jj = 1:N_iter
while diff_ub_lb  > tol_biseciton
    T_val = (T_ub_k+T_lb_k)/2;
    miqp_set_up_problem
    result = gurobi(model, params);

    CPU_time_gurobi = CPU_time_gurobi + result.runtime;
    if isequal(result.status,'OPTIMAL')
        T_ub_k = T_val;
        w_opt_gurobi = result.x;
        T_opt = T_val;
        sucess  = sucess+1;
    else
        T_lb_k = T_val;
    end
    iter = iter+1;
    ub_k = [ub_k;T_ub_k];
    lb_k = [lb_k;T_lb_k];
    T_val_k = [T_val_k;T_opt];
    diff_ub_lb = ub_k(end)-lb_k(end);
    if verbose_gurobi
        fprintf(['---------------------------------------------------------------\n']);
        fprintf(['iter %d, CPU iter: %2.3f, gurobi status: ' result.status ' \n'],iter,result.runtime);
        fprintf(['T_opt= %2.4f lb_k = %2.4f, ub_k = %2.4f.\n'],T_opt,lb_k(end),ub_k(end));
        fprintf(['acc= %2.2e .\n'],abs(lb_k(end)-ub_k(end)));
        fprintf(['---------------------------------------------------------------\n']);
    end
end

T_final_opt = T_val_k(end);
u_opt = w_opt_gurobi(ind_u);
u_opt = reshape(u_opt,n_u,length(u_opt)/n_u);
x_opt = w_opt_gurobi(ind_x);
x_opt = reshape(x_opt,n_x,length(x_opt)/n_x);
y_opt = w_opt_gurobi(ind_y);
y_opt = reshape(y_opt,n_y,length(y_opt)/n_y);

%% stats of solver
if plot_gurobi 
    figure
    subplot(131)
    stairs(lb_k)
    hold on
    stairs(ub_k)
    grid on
    xlabel('iter','Interpreter','latex')
    subplot(132)
    plot(T_val_k,'k')
    xlabel('iter','Interpreter','latex')
    ylabel('$T^*$','Interpreter','latex')
    grid on
    subplot(133)
    semilogy(abs(T_val_k-T_val_k(end)),'k')
    xlabel('iter','Interpreter','latex')
    ylabel('$|T^*-T_k|$','Interpreter','latex')
    grid on
end
