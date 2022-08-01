%% Plot robot kinematics
plot_head = 0;
scale_up = 0;
FontSize = 12;
q0 = [0;0.3;-pi/4;pi/2];
v0 = zeros(4,1);
qx = q0(1);
qz = q0(2);
diff_states = [q0;v0];
p_head_res = full(p_head_fun(diff_states));
p_foot_res = full(p_foot_fun(diff_states));
p_knee_res = full(p_knee_fun(diff_states));
v_foot_res = full(v_foot_fun(diff_states));
v_knee_res = full(v_knee_fun(diff_states));
xx = -0.25:0.5:0.25;
subplot(1,7,1)
hold on
xlabel('$x$','interpreter','latex');
ylabel('$z$','interpreter','latex');

text([qx(1)+0.02],[qz(1)+0.02],'$(q_x,q_z)$','FontSize',FontSize ,'Interpreter','latex')
% thigh
hold on
plot([qx(1),p_knee_res(1,1)],[qz(1),p_knee_res(2,1)],'k','linewidth',6+scale_up);
hold on
% mark angles
direction = -0.4*[qx(1)-p_knee_res(1,1);qz(1)-p_knee_res(2,1)];
plot([qx(1),p_knee_res(1,1)+direction(1)],[qz(1),p_knee_res(2,1)+direction(2)],'k--','linewidth',2);
text(p_knee_res(1)-0.02,p_knee_res(2)-0.07,'$\phi_{\mathrm{knee}}$','FontSize',FontSize ,'Interpreter','latex')

plot([qx(1),qx(1)],[qz(1),qz(1)-0.15],'k--','linewidth',2);
text(qx(1)+0.02,qz(1)-0.1,'$\phi_{\mathrm{hip}}$','FontSize',FontSize ,'Interpreter','latex')
hold on
% shank
plot([p_knee_res(1,1),p_foot_res(1,1)],[p_knee_res(2,1),p_foot_res(2,1)],'k','linewidth',3+scale_up);
% knee
plot([p_knee_res(1,1)],[p_knee_res(2,1)],'ro','markersize',4+scale_up,'markerfacecolor','r');
% foot
% head
hold on
% head - base
plot([qx(1),qx(1)],[qz(1)-0.01,qz(1)+0.03],'k','linewidth',14+scale_up);
% base
plot([qx(1)],[qz(1)],'ro','markersize',9+scale_up,'markerfacecolor','r');
area(xx,0*xx-1,'FaceColor','k','FaceAlpha',0.2);
axis equal
% grid on
xlim([-0.08 0.24]);
ylim([-0.06 0.36]);
% xticklabels([])
% yticklabels([])
