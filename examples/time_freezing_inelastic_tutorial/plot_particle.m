function plot_particle(results, friction)
    %% Extract results
    qx = results.x(1,:);
    qy = results.x(2,:);
    vx = results.x(3,:);
    vy = results.x(4,:);
    t_opt = results.x(5,:);

    %% States
    figure
    subplot(121)
    plot(qx,qy);
    axis equal
    grid on
    xlabel('$q_x$','interpreter','latex');
    ylabel('$q_y$','interpreter','latex');
    subplot(122)
    plot(t_opt,vy);
    hold on
    plot(t_opt,vx);
    grid on
    xlabel('$t$','interpreter','latex');
    ylabel('$v$','interpreter','latex');

    %% plot alpha
    if friction
        alpha1 = results.alpha(1,:);
        alpha2 = results.alpha(2,:);
        alpha3 = results.alpha(3,:);
        theta1 = alpha1+(1-alpha1).*(alpha2);
        alpha_aux = (1-alpha1).*(1-alpha2);
        theta2 = alpha_aux.*(1-alpha3);
        theta3 = alpha_aux.*(alpha3);
        figure;
        subplot(131)
        plot(results.t_grid,[theta1,nan])
        xlabel('$\tau$','interpreter','latex');
        ylabel(['$\theta_1$'],'interpreter','latex');
        grid on
        ylim([-0.1 1.1]);
        subplot(132)
        plot(results.t_grid,[theta2,nan])
        xlabel('$\tau$','interpreter','latex');
        ylabel(['$\theta_2$'],'interpreter','latex');
        grid on
        ylim([-0.1 1.1]);
        subplot(133)
        plot(results.t_grid,[theta3,nan])
        xlabel('$\tau$','interpreter','latex');
        ylabel(['$\theta_3$'],'interpreter','latex');
        grid on
        ylim([-0.1 1.1]);
    else
        alpha1 = results.alpha(1,:);
        alpha2 = results.alpha(2,:);
        theta1 = alpha1+(1-alpha1).*(alpha2);
        theta2 = (1-alpha1).*(1-alpha2);
        figure;
        subplot(121)
        plot(results.t_grid,[theta1,nan])
        xlabel('$\tau$','interpreter','latex');
        ylabel(['$\theta_1$'],'interpreter','latex');
        grid on
        ylim([-0.1 1.1]);
        subplot(122)
        plot(results.t_grid,[theta2,nan])
        xlabel('$\tau$','interpreter','latex');
        ylabel(['$\theta_2$'],'interpreter','latex');
        grid on
        ylim([-0.1 1.1]);
    end
    %% speed of time
    figure
    subplot(121)
    plot(results.t_grid,t_opt)
    hold on
    plot(results.t_grid,results.t_grid,'k--')
    grid on
    xlabel('$\tau$','interpreter','latex');
    ylabel('$t$','interpreter','latex');
    subplot(122)
    stairs(results.s_sot)
    grid on
    xlabel('simulation step','interpreter','latex');
    ylabel('$s$','interpreter','latex');
end
