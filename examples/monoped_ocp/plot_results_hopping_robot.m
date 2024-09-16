import casadi.*
close all;
%%
LineWidth = 3;
make_animation = 0;
plot_kinematics = 0;
label_fontsize = 24;
%% visualization
h = 0.1;
x_res = ocp_solver.get('x');
u_opt = ocp_solver.get('u');
qx = x_res(1,:);
qz = x_res(2,:);
q_phi_hip = x_res(3,:);
q_phi_knee = x_res(4,:);
vx = x_res(5,:);
vz = x_res(6,:);
omega_hip = x_res(7,:);
omega_knee = x_res(8,:);
t_opt = x_res(9,:);
if ~exist('hole_constraint','var')
    hole_constraint = 0;
end
%%
q_res = [qx;qz;q_phi_hip;q_phi_knee];
v_res = [vx;vz;omega_hip;omega_knee];
diff_states = [q_res;v_res];
p_head_res = [];
p_foot_res = [];
p_knee_res = [];
v_foot_res = [];
v_knee_res = [];
% Compute kinematics
for ii = 1:length(diff_states)
    p_head_res = [p_head_res,full(p_head_fun(diff_states(:,ii)))];
    p_foot_res = [p_foot_res,full(p_foot_fun(diff_states(:,ii)))];
    p_knee_res = [p_knee_res,full(p_knee_fun(diff_states(:,ii)))];
    v_foot_res = [v_foot_res,full(v_foot_fun(diff_states(:,ii)))];
    v_knee_res = [v_knee_res,full(v_knee_fun(diff_states(:,ii)))];
end

%% Make gif animation
ymin = min(p_foot_res(2,:))-0.1;
ymax = max(qz)+0.1;
xmin = min(qx)-0.5;
xmax = max(qx)+0.5;
xx = xmin:0.1:xmax;
%
save_gif = 1;
if make_animation
    figure
    for ii = 1:length(qz)
        hold on
        % head - base
        plot([p_head_res(1,ii),qx(ii)],[p_head_res(2,ii),qz(ii)],'k','linewidth',13);
        % base
        plot([qx(ii)],[qz(ii)],'ro','markersize',7,'markerfacecolor','r');
        hold on
        % thigh
        plot([qx(ii),p_knee_res(1,ii)],[qz(ii),p_knee_res(2,ii)],'k','linewidth',5);
        % shank
        plot([p_knee_res(1,ii),p_foot_res(1,ii)],[p_knee_res(2,ii),p_foot_res(2,ii)],'k','linewidth',3);
        % knee
        plot([p_knee_res(1,ii)],[p_knee_res(2,ii)],'bo','markersize',4,'markerfacecolor','b');
        % foot
        plot([p_foot_res(1,ii)],[p_foot_res(2,ii)]+0.005,'ro','markersize',3,'markerfacecolor','r');
        area(xx,0*xx-1,'FaceColor','k','FaceAlpha',0.2);
        if hole_constraint
            for kk= 1:n_holes
                x_nodes = [xc_vec(kk)-a_vec(kk)  xc_vec(kk)+a_vec(kk) xc_vec(kk)+a_vec(kk) xc_vec(kk)-a_vec(kk)];
                y_nodes = [-5 -5 1e-2 1e-2];
                patch(x_nodes,y_nodes,'white','EdgeColor','white');
            end
        end
        axis equal
        xlim([xmin xmax]);
        ylim([ymin ymax]);
        drawnow
        % Capture the plot as an image
        frame = getframe();
        im = frame2im(frame);
        if save_gif
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            if ii == 1
                imwrite(imind,cm,[filename '.gif'],'gif','DelayTime',h,'Loopcount',inf);
            else
                imwrite(imind,cm,[filename '.gif'],'gif','DelayTime',h,'WriteMode','append');
            end
        end
        clf
    end
end
%% Severla frames of the solution
frame_skip = 4;
figure('Renderer', 'painters', 'Position', [100 100 1200 500])
if plot_kinematics
    subplot(1,7,[2:7])
end
area(xx,0*xx-1,'FaceColor','k','FaceAlpha',0.2);
hold on
if hole_constraint
    for kk= 1:n_holes
        x_nodes = [xc_vec(kk)-a_vec(kk)  xc_vec(kk)+a_vec(kk) xc_vec(kk)+a_vec(kk) xc_vec(kk)-a_vec(kk)];
        y_nodes = [-5 -5 0 0];
        patch(x_nodes,y_nodes,'white','EdgeColor','white');
    end
end
hold on
%xlabel('$x$','interpreter','latex');
%ylabel('$z$','interpreter','latex');
plot_head = 1;
for ii = 1:frame_skip:length(qz)
    % thigh
    plot([qx(ii),p_knee_res(1,ii)],[qz(ii),p_knee_res(2,ii)],'k','linewidth',8);
    hold on
    % shank
    plot([p_knee_res(1,ii),p_foot_res(1,ii)],[p_knee_res(2,ii),p_foot_res(2,ii)],'k','linewidth',4);
    % knee
    plot([p_knee_res(1,ii)],[p_knee_res(2,ii)],'ro','markersize',5,'markerfacecolor','r');
    % foot
    hold on
    % head - base
    plot([p_head_res(1,ii),qx(ii)],[p_head_res(2,ii),qz(ii)],'k','linewidth',16);
    % base
    plot([qx(ii)],[qz(ii)],'ro','markersize',10,'markerfacecolor','r');
end

axis equal
%axis tight
xlim([xmin+0.2 xmax-0.2]);
ylim([ymin ymax]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'YTick',[]);
set(gca,'XTick',[]);
%set(gca,'LooseInset',get(gca,'TightInset'));
% static plot of kinemtics
if plot_kinematics
    plot_robot_kinematics
end
if save_figure
    exportgraphics(gca, [filename, '_key_frames.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
end

%% states
LineWidth = LineWidth-0.5;
figure
subplot(221);
plot(t_opt,qx,'Color',0.5*ones(3,1),'LineWidth',LineWidth);
hold on
plot(t_opt,qz,'k','LineWidth',LineWidth);
xlabel('$t$','interpreter','latex', 'Fontsize', label_fontsize);
ylabel('Position of the base','interpreter','latex', 'Fontsize', label_fontsize);
legend({'$x_{\textrm{base}}$','$z_{\textrm{base}}$'},'interpreter','latex','Location','best');
grid on
subplot(223);
plot(t_opt,q_phi_hip,'k','LineWidth',LineWidth);
hold on
plot(t_opt,q_phi_knee,'Color',0.5*ones(3,1),'LineWidth',LineWidth);
ylim([[min([q_phi_hip,q_phi_knee])]-1  [max([q_phi_hip,q_phi_knee])]+1  ])
% plot(tt,0*q_knee(ind)-pi/2*1.05,'r');
% plot(tt,0*q_knee(ind)+pi/2*1.05,'r');
xlabel('$t$','interpreter','latex', 'Fontsize', label_fontsize);
ylabel('$\phi$','interpreter','latex', 'Fontsize', label_fontsize);
legend({'$\phi_{\textrm{hip}}$','$\phi_{\textrm{knee}}$'},'interpreter','latex','Location','best');
grid on
subplot(222);
plot(t_opt,vx,'Color',0.5*ones(3,1),'LineWidth',LineWidth);
hold on
plot(t_opt,vz,'k','LineWidth',LineWidth);
xlabel('$t$','interpreter','latex', 'Fontsize', label_fontsize);
ylabel('Velocity of the base','interpreter','latex', 'Fontsize', label_fontsize);
legend({'$v^x_{\textrm{base}}$','$v^z_{\textrm{base}}$'},'interpreter','latex','Location','best');
grid on
subplot(224);
plot(t_opt,omega_hip,'k','LineWidth',LineWidth);
hold on
plot(t_opt,omega_knee,'Color',0.5*ones(3,1),'LineWidth',LineWidth);
ylim([[min([omega_knee,omega_hip])]-1  [max([omega_knee,omega_hip])]+1  ])
xlabel('$t$','interpreter','latex', 'Fontsize', label_fontsize);
ylabel('$\omega$','interpreter','latex', 'Fontsize', label_fontsize);
legend({'$\omega_{\textrm{hip}}$','$\omega_{\textrm{knee}}$'},'interpreter','latex','Location','best');
grid on

%% knee and foot
figure
subplot(221);
plot(t_opt,p_knee_res(1,:),'Color',0.5*ones(3,1),'LineWidth',LineWidth);
hold on
plot(t_opt,p_foot_res(1,:),'k','LineWidth',LineWidth);
xlabel('$t$','interpreter','latex', 'Fontsize', label_fontsize);
ylabel('Horizontal Positions','interpreter','latex', 'Fontsize', label_fontsize);
legend({'$p^x_{\mathrm{knee}}$','$p^x_{\mathrm{foot}}$'},'interpreter','latex','location','best');
grid on

subplot(222);
plot(t_opt,p_knee_res(2,:),'Color',0.5*ones(3,1),'LineWidth',LineWidth);
hold on
plot(t_opt,p_foot_res(2,:),'k','LineWidth',LineWidth);
xlabel('$t$','interpreter','latex', 'Fontsize', label_fontsize);
ylabel('Vertical Positions','interpreter','latex', 'Fontsize', label_fontsize);
legend({'$p^y_{\mathrm{knee}}$','$p^y_{\mathrm{foot}}$'},'interpreter','latex','location','best');
grid on

subplot(223);
plot(t_opt,v_knee_res(1,:),'Color',0.5*ones(3,1),'LineWidth',LineWidth);
hold on
plot(t_opt,v_foot_res(1,:),'k','LineWidth',LineWidth);
xlabel('$t$','interpreter','latex', 'Fontsize', label_fontsize);
ylabel('Horizontal Velocities','interpreter','latex', 'Fontsize', label_fontsize);
legend({'$v_{\mathrm{knee}}$','$v_{\mathrm{foot}}$'},'interpreter','latex','location','best');
grid on

subplot(224);
plot(t_opt,v_knee_res(2,:),'Color',0.5*ones(3,1),'LineWidth',LineWidth);
hold on
plot(t_opt,v_foot_res(2,:),'k','LineWidth',LineWidth);
xlabel('$t$','interpreter','latex', 'Fontsize', label_fontsize);
ylabel('Vertical Velocities','interpreter','latex', 'Fontsize', label_fontsize);
legend({'$v^y_{\mathrm{knee}}$','$v^y_{\mathrm{foot}}$'},'interpreter','latex','location','best');
grid on

if save_figure
    exportgraphics(gca, [filename, '_kinematics.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
end
LineWidth = LineWidth+0.5;
%% Plot controls
if model.dims.n_u > 0
    % refine control;
    figure('Renderer', 'painters')
    t_grid_u = ocp_solver.get_control_grid();
    stairs(t_grid_u,[u_opt(1,:),nan],"Color",'k','LineWidth',LineWidth);
    hold on
    stairs(t_grid_u,[u_opt(2,:),nan],"Color",0.5*ones(3,1),'LineWidth',LineWidth);
    xlabel('$t$','interpreter','latex', 'Fontsize', label_fontsize);
    ylabel('$u(t)$','interpreter','latex', 'Fontsize', label_fontsize);
    grid on
    legend({'$u_{\mathrm{hip}}(t)$','$u_{\mathrm{knee}}(t)$'},'interpreter','latex','location','best');
    xlim([0 t_grid_u(end)+1e-6]);
    if save_figure
        exportgraphics(gca, [filename, '_controls.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
    end

    try
        figure
        s_sot = ocp_solver.get('sot');
        stairs(t_grid_u,[s_sot,nan],"Color",0.2*ones(3,1),'LineWidth',LineWidth);
        xlabel('$t$','interpreter','latex', 'Fontsize', label_fontsize);
        ylabel('$s(t)$','interpreter','latex', 'Fontsize', label_fontsize);
        xlim([0 t_grid_u(end)]);
        grid on
        if save_figure
            exportgraphics(gca, [filename, '_sot.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
        end

    catch
    end
end


%% constraint drif, normal and tangential velocity
f_c = model.f_c;
tangent1 = model.J_tangent;
nabla_q_f_c = model.J_normal;
f_c_fun = Function('f_c_fun',{q},{f_c});
nabla_q_f_c_fun = Function('nabla_q_f_c_fun',{q},{nabla_q_f_c});
tangent1_fun = Function('tangent1_fun',{q},{tangent1});
constraint_drift = [];
for ii = 1:length(q_res)
    temp = [full(nabla_q_f_c_fun(q_res(:,ii)))'*v_res(:,ii);...
        full(f_c_fun(q_res(:,ii)));...
        full(tangent1_fun(q_res(:,ii)))'*v_res(:,ii)];
    constraint_drift = [constraint_drift, temp];
end
T = x_res(end,end);

figure
subplot(311)
plot(t_opt,constraint_drift(2,:),'k','LineWidth',LineWidth);
ylim([-0.1 max(constraint_drift(2,:))+0.1]);
xlim([0 T]);

grid on
xlabel('$t$','interpreter','latex', 'Fontsize', label_fontsize);
ylabel('$f_c(q)$','interpreter','latex', 'Fontsize', label_fontsize);
subplot(312)
plot(t_opt,constraint_drift(1,:),'k','LineWidth',LineWidth);
grid on
xlabel('$t$','interpreter','latex', 'Fontsize', label_fontsize);
ylabel('$n(q)^\top v$','interpreter','latex', 'Fontsize', label_fontsize);
xlim([0 T]);
ylim([ min(constraint_drift(1,:))-0.1 max(constraint_drift(1,:))+0.1]);
subplot(313)
plot(t_opt,constraint_drift(3,:),'k','LineWidth',LineWidth);
grid on
xlabel('$t$','interpreter','latex', 'Fontsize', label_fontsize);
ylabel('$t_1(q)^\top v$','interpreter','latex', 'Fontsize', label_fontsize);
ylim([min(constraint_drift(3,:))-0.1 max(constraint_drift(3,:))+0.1]);
xlim([0 T]);
if save_figure
    exportgraphics(gca, [filename, '_drifts.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
end

