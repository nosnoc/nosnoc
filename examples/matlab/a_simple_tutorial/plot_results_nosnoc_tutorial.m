%%  Plot results
t_grid = results.t_grid;
t_grid_u = results.t_grid_u;
q_opt=results.x_opt(1,:);
v_opt=results.x_opt(2,:);
u_opt=results.u_opt;

if settings.time_optimal_problem
    fprintf('Final time: %2.4f s.\n',results.T_opt)
else
    fprintf('Objective value time: %2.4f s.\n',results.f_opt)
end

subplot(311)
plot(t_grid,q_opt)
xlabel('$t$','Interpreter','latex')
xlabel('$q(t)$','Interpreter','latex')
grid on
subplot(312)
plot(t_grid,v_opt)
yline(v_max,'r--')
yline(v_trash_hold,'k--')
ylim(1.2*[0 v_max]);
xlabel('$t$','Interpreter','latex')
xlabel('$v(t)$','Interpreter','latex')
grid on
subplot(313)
stairs(t_grid_u,[u_opt,nan])
xlabel('$t$','Interpreter','latex')
xlabel('$u(t)$','Interpreter','latex')
yline(-u_max,'r--')
yline(u_max,'r--')
ylim(1.2*[-u_max u_max]);
grid on
%%
% T_star=13+1/3;
% error = norm(results.T_opt-T_star);
% fprintf('Numerical error %2.2e.\n',error)