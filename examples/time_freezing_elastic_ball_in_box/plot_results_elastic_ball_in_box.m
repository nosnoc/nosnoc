file_name1 = 'binbx_clock_state1';
file_name2 = 'binbx_state_and_controls_1';
file_name3 = 'binbx_trajectory1';

%  Colors
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
organe = [0.9290 0.6940 0.1250];
grey = [0.85 0.85 0.85];

qx_opt = x_opt(1,:);
qy_opt = x_opt(2,:);
vx_opt = x_opt(3,:);
vy_opt = x_opt(4,:);
t_opt = x_opt(5,:);
%%

figure('Renderer', 'painters', 'Position', [100 100 800 200])
subplot(121)
plot(t_grid,t_opt,'k-','linewidth',1.5)
xlim([0 t_opt(end)])
hold on
t_fin = t_opt(end);
plot(0:0.1:t_fin,0:0.1:t_fin,'r-')
xlabel('Numerical time - $\tau$','Interpreter','latex');
ylabel('Physical time - $t$','Interpreter','latex');
hold on
for ii= 1:N_finite_elements(1):length(t_grid)
    xline(t_grid(ii),'k--')
end
s_sot_opt = results.w_opt(model.ind_sot);
subplot(122)
stairs(t_grid_u,[nan;s_sot_opt],'k','linewidth',1.5);
ylim([0.0 max(s_sot_opt)+1])
xlabel('$\tau$','Interpreter','latex');
ylabel('$s(\tau)$','Interpreter','latex');
grid on
xlim([0 t_opt(end)])
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
eval(['print -dpdf -painters ' file_name1])

%% 
linewidth = 1.5;
ind_t = diff(t_opt)/h>0.2;
ind_t = [1,theta_opt(1,:)]>0.2;
figure('Renderer', 'painters', 'Position', [100 100 700 850])
subplot(321)
plot(t_grid,qx_opt,'linewidth',linewidth)
hold on
plot(t_grid,qy_opt,'linewidth',linewidth)
plot(t_grid,R*sin(omega*(t_opt)+alpha0),':','Color',blue,'linewidth',1.2);
plot(t_grid,R*cos(omega*(t_opt)+alpha0),':','Color',red,'linewidth',1.2);
grid on
xlabel('$\tau$','Interpreter','latex');
ylabel('$q(\tau)$','Interpreter','latex');
xlim([0 t_opt(end)])
ylim ([-2 3])
l = legend({'$q_1(\tau)$','$q_2(\tau)$','$q^{\mathrm{ref}}_1(\tau)$','$q^{\mathrm{ref}}_2(\tau)$'},'Interpreter','latex','location','northeast');
l.NumColumns = 2;

subplot(323)
plot(t_grid,vx_opt,'linewidth',linewidth)
hold on
grid on
plot(t_grid,vy_opt,'linewidth',linewidth)
xlabel('$\tau$','Interpreter','latex');
ylabel('$v(\tau)$','Interpreter','latex');
legend({'$v_1(\tau)$','$v_2(\tau)$'},'Interpreter','latex','location','north');
xlim([0 t_opt(end)])
ylim ([-11 11])

subplot(325)
stairs(t_grid_u,[nan,u_opt(1,:)],'linewidth',linewidth)
hold on
grid on
stairs(t_grid_u,[nan,u_opt(2,:)],'linewidth',linewidth)
xlabel('$\tau$','Interpreter','latex');
ylabel('$u(\tau)$','Interpreter','latex');
legend({'$u_1(\tau)$','$u_2(\tau)$'},'Interpreter','latex','location','north');
xlim([0 t_opt(end)])
ylim([-50 60])


subplot(322)
plot(t_opt(ind_t),qx_opt(ind_t),'linewidth',linewidth)
hold on
plot(t_opt(ind_t),qy_opt(ind_t),'linewidth',linewidth)
plot(t_opt(ind_t),R*sin(omega*(t_opt(ind_t))+alpha0),':','Color',blue,'linewidth',1.2);
plot(t_opt(ind_t),R*cos(omega*(t_opt(ind_t))+alpha0),':','Color',red,'linewidth',1.2);
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$q(t$','Interpreter','latex');
l = legend({'$q_1(t)$','$q_2(t)$','$q^{\mathrm{ref}}_1(t)$','$q^{\mathrm{ref}}_2(t)$'},'Interpreter','latex','location','northeast');
l.NumColumns = 2;
xlim([0 t_opt(end)])
ylim ([-2 3])


subplot(324)
plot(t_opt,vx_opt,'linewidth',linewidth)
hold on
grid on
plot(t_opt,vy_opt,'linewidth',linewidth)
xlabel('$t$','Interpreter','latex');
ylabel('$v(t)$','Interpreter','latex');
legend({'$v_1(t)$','$v_2(t)$'},'Interpreter','latex','location','northeast');
xlim([0 t_opt(end)])
ylim ([-11 11])

subplot(326)
stairs(t_opt(1:N_finite_elements(1):end),[nan,u_opt(1,:)],'linewidth',linewidth)
hold on
grid on
stairs(t_opt(1:N_finite_elements(1):end),[nan,u_opt(2,:)],'linewidth',linewidth)
xlabel('$t$','Interpreter','latex');
ylabel('$u(t)$','Interpreter','latex');
legend({'$u_1(t)$','$u_2(t)$'},'Interpreter','latex','location','north');
xlim([0 t_opt(end)])
ylim([-50 60])

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
eval(['print -dpdf -painters ' file_name2])
%%
figure('Renderer', 'painters', 'Position', [100 100 350 300])
plot(qx_opt(ind_t),qy_opt(ind_t),'LineWidth',2);
hold on
fimplicit(@(x,y) (x-qx_c).^2+(y-qy_c).^2-R^2, [-2 2],'Color',blue,'LineStyle',':','LineWidth',1.5);
plot(qx_opt(1),qy_opt(1),'r.','MarkerSize',12);
plot(qx_opt(end),qy_opt(end),'rx','MarkerSize',7);
try
    plot([a_left, a_right],[b_bottom, b_bottom],'Color','k','LineWidth',1.5);
    plot([a_left, a_right],[b_top, b_top],'Color','k','LineWidth',1.5);
    plot([a_left, a_left],[b_bottom, b_top],'Color','k','LineWidth',1.5);
    plot([a_right, a_right],[b_bottom, b_top],'Color','k','LineWidth',1.5);
catch
    yline(0,'r');
end
grid on
legend({'$q$','$q^{\mathrm{ref}}$','Starpoint','Endpoint'},'interpreter','latex','Location','best');
xlim([-1.2 1.2])
ylim([-1.2 1.2])
xlabel('$q_1$','interpreter','latex');
ylabel('$q_2$','interpreter','latex');
axis equal

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
eval(['print -dpdf -painters ' file_name3])