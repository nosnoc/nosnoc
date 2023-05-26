%%
figure
subplot(211)
plot(t_grid,x(1,:));
hold on
plot(t_grid,x(2,:));
plot(t_grid,x(3,:));
xlabel('$t$','interpreter','latex');
ylabel('$q(t)$','interpreter','latex');
legend({'$q_1(t)$','$q_2(t)$','$q_3(t)$'},'interpreter','latex');
grid on
subplot(212)
plot(t_grid,x(4,:));
hold on
plot(t_grid,x(5,:));
plot(t_grid,x(6,:));
xlabel('$t$','interpreter','latex');
ylabel('$v(t)$','interpreter','latex');
legend({'$v_1(t)$','$v_2(t)$','$v_3(t)$'},'interpreter','latex');
grid on

switch dcs_mode
    case 'Stewart'
        figure
        subplot(311)
        plot(t_grid,[theta(1,:),nan])
        hold on
        grid on
        plot(t_grid,[theta(2,:),nan])
        xlabel('$t$','interpreter','latex');
        ylabel('$\theta(t)$','interpreter','latex');
        legend({'$\theta_1(t)$','$\theta_2(t)$'},'interpreter','latex');
        subplot(312)
        plot(t_grid,[theta(3,:),nan])
        hold on
        grid on
        plot(t_grid,[theta(4,:),nan])
        xlabel('$t$','interpreter','latex');
        ylabel('$\theta(t)$','interpreter','latex');
        legend({'$\theta_3(t)$','$\theta_4(t)$'},'interpreter','latex');
        subplot(313)
        plot(t_grid,[theta(5,:),nan])
        hold on
        grid on
        plot(t_grid,[theta(6,:),nan])
        xlabel('$t$','interpreter','latex');
        ylabel('$\theta(t)$','interpreter','latex');
        legend({'$\theta_5(t)$','$\theta_6(t)$'},'interpreter','latex');
    case 'Step'
        figure
        subplot(311)
        plot(t_grid,[alpha(1,:),nan])
        grid on
        xlabel('$t$','interpreter','latex');
        ylabel('$\alpha_1(t)$','interpreter','latex');
        subplot(312)
        plot(t_grid,[alpha(2,:),nan])
        grid on
        xlabel('$t$','interpreter','latex');
        ylabel('$\alpha_2(t)$','interpreter','latex');
        subplot(313)
        plot(t_grid,[alpha(3,:),nan])
        grid on
        xlabel('$t$','interpreter','latex');
        ylabel('$\alpha_3(t)$','interpreter','latex');

end

%%
figure
switch dcs_mode
    case 'Stewart'
        subplot(311)
        plot(t_grid,[lam(1,:),nan])
        hold on
        grid on
        plot(t_grid,[lam(2,:),nan])
        xlabel('$t$','interpreter','latex');
        ylabel('$\lambda(t)$','interpreter','latex');
        legend({'$\lambda_1(t)$','$\lambda_2(t)$'},'interpreter','latex');
        subplot(312)
        plot(t_grid,[lam(3,:),nan])
        hold on
        grid on
        plot(t_grid,[lam(4,:),nan])
        xlabel('$t$','interpreter','latex');
        ylabel('$\lambda(t)$','interpreter','latex');
        legend({'$\lambda_3(t)$','$\lambda_4(t)$'},'interpreter','latex');
        subplot(313)
        plot(t_grid,[lam(5,:),nan])
        hold on
        grid on
        plot(t_grid,[lam(6,:),nan])
        xlabel('$t$','interpreter','latex');
        ylabel('$\lambda(t)$','interpreter','latex');
        legend({'$\lambda_5(t)$','$\lambda_6(t)$'},'interpreter','latex');
        %
        figure
        plot(t_grid,[mu(1,:),nan])
        hold on
        grid on
        plot(t_grid,[mu(2,:),nan])
        plot(t_grid,[mu(3,:),nan])
        xlabel('$t$','interpreter','latex');
        ylabel('$\mu(t)$','interpreter','latex');
        legend({'$\mu_1(t)$','$\mu_2(t)$','$\mu_3(t)$'},'interpreter','latex');
    case 'Step'
        subplot(311)
        plot(t_grid,[lambda_n(1,:),nan])
        hold on
        grid on
        plot(t_grid,[lambda_p(1,:),nan])
        xlabel('$t$','interpreter','latex');
        ylabel('$\lambda_1(t)$','interpreter','latex');
        legend({'$\lambda_{0,1}(t)$','$\lambda_{1,1}(t)$'},'interpreter','latex');
        subplot(312)
        plot(t_grid,[lambda_n(2,:),nan])
        hold on
        grid on
        plot(t_grid,[lambda_p(2,:),nan])
        xlabel('$t$','interpreter','latex');
        ylabel('$\lambda_2(t)$','interpreter','latex');
        legend({'$\lambda_{0,2}(t)$','$\lambda_{1,2}(t)$'},'interpreter','latex');
        subplot(313)
        plot(t_grid,[lambda_n(3,:),nan])
        hold on
        grid on
        plot(t_grid,[lambda_p(3,:),nan])
        xlabel('$t$','interpreter','latex');
        ylabel('$\lambda_3(t)$','interpreter','latex');
        legend({'$\lambda_{0,3}(t)$','$\lambda_{1,3}(t)$'},'interpreter','latex');
end
%%
figure
stairs(results.h,'k')
xlabel('finite element','interpreter','latex');
ylabel('$h_{n}$','interpreter','latex');
