%%
figure
subplot(211)
plot(t_grid,x_res(1,:));
hold on
plot(t_grid,x_res(2,:));
plot(t_grid,x_res(3,:));
xlabel('$t$','interpreter','latex');
ylabel('$q(t)$','interpreter','latex');
legend({'$q_1(t)$','$q_2(t)$','$q_3(t)$'},'interpreter','latex');
grid on
subplot(212)
plot(t_grid,x_res(4,:));
hold on
plot(t_grid,x_res(5,:));
plot(t_grid,x_res(6,:));
xlabel('$t$','interpreter','latex');
ylabel('$v(t)$','interpreter','latex');
legend({'$v_1(t)$','$v_2(t)$','$v_3(t)$'},'interpreter','latex');
grid on

figure
subplot(311)
plot(t_grid,[theta_res(1,:),nan])
hold on
grid on
plot(t_grid,[theta_res(2,:),nan])
xlabel('$t$','interpreter','latex');
ylabel('$\theta(t)$','interpreter','latex');
legend({'$\theta_1(t)$','$\theta_2(t)$'},'interpreter','latex');
subplot(312)
plot(t_grid,[theta_res(3,:),nan])
hold on
grid on
plot(t_grid,[theta_res(4,:),nan])
xlabel('$t$','interpreter','latex');
ylabel('$\theta(t)$','interpreter','latex');
legend({'$\theta_3(t)$','$\theta_4(t)$'},'interpreter','latex');
subplot(313)
plot(t_grid,[theta_res(5,:),nan])
hold on
grid on
plot(t_grid,[theta_res(6,:),nan])
xlabel('$t$','interpreter','latex');
ylabel('$\theta(t)$','interpreter','latex');
legend({'$\theta_5(t)$','$\theta_6(t)$'},'interpreter','latex');

%%
figure
subplot(311)
plot(t_grid,[lambda_res(1,:),nan])
hold on
grid on
plot(t_grid,[lambda_res(2,:),nan])
xlabel('$t$','interpreter','latex');
ylabel('$\lambda(t)$','interpreter','latex');
legend({'$\lambda_1(t)$','$\lambda_2(t)$'},'interpreter','latex');
subplot(312)
plot(t_grid,[lambda_res(3,:),nan])
hold on
grid on
plot(t_grid,[lambda_res(4,:),nan])
xlabel('$t$','interpreter','latex');
ylabel('$\lambda(t)$','interpreter','latex');
legend({'$\lambda_3(t)$','$\lambda_4(t)$'},'interpreter','latex');
subplot(313)
plot(t_grid,[lambda_res(5,:),nan])
hold on
grid on
plot(t_grid,[lambda_res(6,:),nan])
xlabel('$t$','interpreter','latex');
ylabel('$\lambda(t)$','interpreter','latex');
legend({'$\lambda_5(t)$','$\lambda_6(t)$'},'interpreter','latex');
%%
figure
plot(t_grid,[mu_res(1,:),nan])
hold on
grid on
plot(t_grid,[mu_res(2,:),nan])
plot(t_grid,[mu_res(3,:),nan])
xlabel('$t$','interpreter','latex');
ylabel('$\mu(t)$','interpreter','latex');
legend({'$\mu_1(t)$','$\mu_2(t)$','$\mu_3(t)$'},'interpreter','latex');
%%
figure
stairs(results.h_vec,'k')
xlabel('finite element','interpreter','latex');
ylabel('$h_{n}$','interpreter','latex');
%%
figure
semilogy(stats.complementarity_stats+1e-20,'k','LineWidth',1.5)
xlabel('integration step n','interpreter','latex');
ylabel('comp residual','interpreter','latex');
grid on


%% Strong stationarity
% figure
% plot(theta1_opt+lambda1_opt)
% hold on
% plot(theta2_opt+lambda2_opt)
% plot(theta3_opt+lambda3_opt)
% grid on
% yline(0,'k')
% ylabel('$\theta+\lambda>0$','interpreter','latex');
% 
% a = min(theta1_opt+lambda1_opt);
% b = min(theta2_opt+lambda2_opt);
% c = min(theta3_opt+lambda3_opt);
% biactivity = min([a;b;c])

