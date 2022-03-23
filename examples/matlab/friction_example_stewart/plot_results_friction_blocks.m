%% Read and plot Result 
% differential states
for i = 1:n_x
    eval( ['x' num2str(i) '_opt = results.x_res(' num2str(i) ',:);']);
end
for i = 1:n_theta
    eval( ['theta' num2str(i) '_opt = results.theta_res(' num2str(i) ',:);']);
end

for i = 1:n_theta
    eval( ['lambda' num2str(i) '_opt = results.lambda_res(' num2str(i) ',:);']);
end

for i = 1:n_simplex
    eval( ['mu' num2str(i) '_opt = results.mu_res(' num2str(i) ',:);']);
end
%%
figure
subplot(211)
plot(t_grid,x1_opt);
hold on
plot(t_grid,x2_opt);
plot(t_grid,x3_opt);
xlabel('$t$','interpreter','latex');
ylabel('$q(t)$','interpreter','latex');
legend({'$q_1(t)$','$q_2(t)$','$q_3(t)$'},'interpreter','latex');
grid on
subplot(212)
plot(t_grid,x4_opt);
hold on
plot(t_grid,x5_opt);
plot(t_grid,x6_opt);
xlabel('$t$','interpreter','latex');
ylabel('$v(t)$','interpreter','latex');
legend({'$v_1(t)$','$v_2(t)$','$v_3(t)$'},'interpreter','latex');
grid on
% figure
% plot(t_grid,x7_opt);
% xlabel('$t$','interpreter','latex');
% ylabel('$t(t)$','interpreter','latex');
% grid on
% N_switches = sum((abs(diff(h_opt)))>1e-4);
% fprintf('Number of Switches is %d. \n',N_switches);
figure
subplot(311)
plot(t_grid,[theta1_opt,nan])
hold on
grid on
plot(t_grid,[theta2_opt,nan])
xlabel('$t$','interpreter','latex');
ylabel('$\theta(t)$','interpreter','latex');
legend({'$\theta_1(t)$','$\theta_2(t)$'},'interpreter','latex');
subplot(312)
plot(t_grid,[theta3_opt,nan])
hold on
grid on
plot(t_grid,[theta4_opt,nan])
xlabel('$t$','interpreter','latex');
ylabel('$\theta(t)$','interpreter','latex');
legend({'$\theta_3(t)$','$\theta_4(t)$'},'interpreter','latex');
subplot(313)
plot(t_grid,[theta5_opt,nan])
hold on
grid on
plot(t_grid,[theta6_opt,nan])
xlabel('$t$','interpreter','latex');
ylabel('$\theta(t)$','interpreter','latex');
legend({'$\theta_5(t)$','$\theta_6(t)$'},'interpreter','latex');

%%
figure
subplot(311)
plot(t_grid,[lambda1_opt,nan])
hold on
grid on
plot(t_grid,[lambda2_opt,nan])
xlabel('$t$','interpreter','latex');
ylabel('$\lambda(t)$','interpreter','latex');
legend({'$\lambda_1(t)$','$\lambda_2(t)$'},'interpreter','latex');
subplot(312)
plot(t_grid,[lambda3_opt,nan])
hold on
grid on
plot(t_grid,[lambda4_opt,nan])
xlabel('$t$','interpreter','latex');
ylabel('$\lambda(t)$','interpreter','latex');
legend({'$\lambda_3(t)$','$\lambda_4(t)$'},'interpreter','latex');
subplot(313)
plot(t_grid,[lambda5_opt,nan])
hold on
grid on
plot(t_grid,[lambda6_opt,nan])
xlabel('$t$','interpreter','latex');
ylabel('$\lambda(t)$','interpreter','latex');
legend({'$\lambda_5(t)$','$\lambda_6(t)$'},'interpreter','latex');
%%
figure
plot(t_grid,[mu1_opt,nan])
hold on
grid on
plot(t_grid,[mu2_opt,nan])
plot(t_grid,[mu3_opt,nan])
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

if model.N_stages == 2
    number_of_switches = sum(abs(diff(h_vec))>1e-6)
end

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

