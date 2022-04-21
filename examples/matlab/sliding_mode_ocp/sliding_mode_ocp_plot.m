
%% controls
unfold_struct(model,'caller');
unfold_struct(results,'caller');
u1_opt = u_opt(1,:);
u2_opt = u_opt(2,:);



x1_opt = x_opt(1,:);
x2_opt = x_opt(2,:);
if linear_control
    v1_opt = x_opt(3,:);
    v2_opt = x_opt(4,:);
end

figure
stairs(t_grid(1:N_finite_elements(1):end),[u1_opt,nan])
hold on
stairs(t_grid(1:N_finite_elements(1):end),[u2_opt,nan])
xlabel('$t$','interpreter','latex');
ylabel('$u(t)$','interpreter','latex');
grid on
xlim([0 model.T])
ylim([-u_max*1.1 u_max*1.1])
%%

if linear_control
    figure
    subplot(211)
    plot(t_grid,v1_opt);
    hold on
    plot(t_grid,v2_opt);
    xlabel('$t$','interpreter','latex');
    ylabel('$v(t)$','interpreter','latex');
    for ii=1:N_stages
        xline(ii*h,'k--')
    end
    grid on
    xlim([0 model.T])
    ylim([-3 3])
    subplot(212)
    stairs(t_grid(1:N_finite_elements(1):end),[u1_opt,nan])
    hold on
    stairs(t_grid(1:N_finite_elements(1):end),[u2_opt,nan])
    xlabel('$t$','interpreter','latex');
    ylabel('$u(t)$','interpreter','latex');
    grid on
    xlim([0 model.T])

end
%% plots
figure
plot(t_grid,x1_opt);
hold on
plot(t_grid,x2_opt);
xlabel('$t$','interpreter','latex');
ylabel('$x(t)$','interpreter','latex');
% for ii=1:N_stages
%     xline(ii*h,'k--')
% end

ctrl_grid  = t_grid(1:N_finite_elements(1):end);
for ii=1:N_stages
    xline(ctrl_grid(ii+1),'k--')
end

for ii=1:length(h_opt)
    xline(t_grid(ii+1),'k:')
end

% grid on
axis equal
legend({'$x_1(t)$','$x_2(t)$','Control Intervals'},'interpreter','latex');
xlim([0 model.T])
%%
figure
plot(x1_opt,x2_opt,'LineWidth',2);
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
grid on
axis equal
if illustrate_regions
    hold on
    p = 2; a = -0.1; a1 = 0; 
    b = -0.05; q = 3;
    t2 = -5:0.01:5;
    plot(-a*(t2-a1).^p,t2,'k')
    hold on
    t1 = t2;
    plot(t1,-b*t1.^q,'k')
    grid on
    axis equal
    xlim([-1.5 1.5])
    ylim([-1.5 1.5])
end



% %%
% if piecewise_equidistant_grid
%     figure
%     nu_opt = full(nu_fun(w_opt))+1e-16;
%     subplot(211)
%     semilogy(tgrid,[nan;(nu_opt);nan]);
%     grid on
%     xlabel('$t$','interpreter','latex');
%     ylabel('$\nu_k(t)$','interpreter','latex');
%     subplot(212)
%     plot(tgrid,[nan;rho_h*tanh(nu_opt*10);nan]);
%     xlabel('$t$','interpreter','latex');
%     ylabel('$\rho_h \tanh(\nu_k(t))$','interpreter','latex');
%     grid on
% end

