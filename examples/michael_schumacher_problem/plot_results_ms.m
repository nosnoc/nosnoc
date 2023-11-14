unfold_struct(model,'caller')
unfold_struct(problem_options,'caller')
%unfold_struct(solver_options,'caller')
unfold_struct(results,'caller')
tgrid = linspace(0, T, N_stages+1);
tgrid_z = linspace(0, T, N_stages);
%% read solutions

u1_opt = results.u(1,:);
u2_opt = results.u(2,:);

x1_opt = results.x(1,:);
x2_opt = results.x(2,:);
x3_opt = results.x(3,:);
x4_opt = results.x(4,:);
x5_opt = results.x(5,:);
%% controls
figure
stairs(t_grid(1:N_finite_elements(1):end),[u1_opt,nan])
hold on
stairs(t_grid(1:N_finite_elements(1):end),[u2_opt,nan])
xlabel('$t$','interpreter','latex');
ylabel('$u(t)$','interpreter','latex');
ylim([-2.2 2.2])
grid on

% %% time plot
% subplot(211)
% plot(t_grid,x1_opt);
% hold on
% plot(t_grid,x2_opt);
% xlabel('$t$','interpreter','latex');
% ylabel('$q$','interpreter','latex');
% grid on
% subplot(212)
% plot(t_grid,x3_opt);
% hold on
% plot(t_grid,x4_opt);
% xlabel('$t$','interpreter','latex');
% ylabel('$v$','interpreter','latex');
% grid on
%% plots
figure
subplot(211)
plot(x1_opt,x2_opt);
hold on
plot(x1_opt,x2_opt,'r.');
xlabel('$q_x$','interpreter','latex');
ylabel('$q_y$','interpreter','latex');
hold on
grid on
plot(q_target(1),q_target(2),'rx','MarkerSize',6)
xx = linspace(0,q_target(1),1e2);

switch path_constraint
    case 'none'
        %         model.g_ineq = [];
    case 'linear'
        plot(xx ,(xx)-track_width,'k');
        plot(xx ,(xx)+track_width,'k');
    case 'nonlinear'
        plot(xx ,sin(omega*xx)-track_width,'k');
        plot(xx ,sin(omega*xx)+track_width,'k');
    case 'track'
        % figure
        arg1 = xx-pi;
        arg2 = xx-2*pi;
        % sig = 1e-3;
        step1 = 0.5*(1+tanh(arg1/sig));
        step2 = 0.5*(1+tanh(arg2/sig));
        yy = sin(xx).*(1-step1)+(pi-xx).*step1.*(1-step2)+(-pi-sin(xx)).*step2;
        plot(xx ,yy-track_width,'k');
        hold on
        plot(xx ,yy+track_width,'k');
    case 'chicane'
        yy = (chicane_width)+chicane_width*tanh(chicane_tightness*(xx-q_target(1)/2));
        plot(xx ,yy-track_width,'k');
        hold on
        plot(xx ,yy+track_width,'k');
end
grid on
axis equal
% xlim([0 3*pi])

subplot(212)
stairs(x1_opt(1:N_finite_elements(1):end),[u1_opt,nan]);
hold on
stairs(x1_opt(1:N_finite_elements(1):end),[u2_opt,nan]);
xlabel('$u(q_x)$','interpreter','latex');
ylabel('$q_x$','interpreter','latex');
legend({'$a(x)$','$s(x)$'},'interpreter','latex');
ylim([-2.2 2.2])
grid on


%%  Normal Force, Veolocity,
v_tangent = [];
v_normal  = [];
v_opt = [x3_opt;x4_opt];
tangent = [cos(x5_opt);sin(x5_opt)];
normal = [-sin(x5_opt);cos(x5_opt)];
for ii = 1 :length(v_opt);
    v_tangent = [v_tangent,tangent(:,ii)'*v_opt(:,ii)];
    v_normal  = [v_normal,normal(:,ii)'*v_opt(:,ii)];
end
figure
 plot(t_grid,v_normal);
hold on
plot(t_grid,v_tangent);
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
legend({'$n^\top v$','$t^\top v$'},'interpreter','latex');
%%
%
% hold on
% %
% %%
% theta2_opt = [theta2_opt(1);theta2_opt];
%
% %% Friction force
% theta1_opt = [theta1_opt(1);theta1_opt];
% figure
% subplot(211)
% plot(x1_opt,Friction_max.*theta1_opt -Friction_max.*theta2_opt);
% xlabel('$q_x$','interpreter','latex');
% ylabel('$F$','interpreter','latex');
% hold on
% plot(x1_opt,normal_velocity,'--');
% grid on
% subplot(212)
% plot(tgrid,Friction_max.*theta1_opt -Friction_max.*theta2_opt);
% xlabel('$t$','interpreter','latex');
% ylabel('$F$','interpreter','latex');
% hold on
% plot(tgrid,normal_velocity,'--');
% grid on
% %%
% figure
% subplot(211)
% plot(tgrid,normal_velocity);
% xlabel('$t$','interpreter','latex');
% ylabel('$v$','interpreter','latex');
% grid on
% subplot(212)
%   for i = m_ind_vec(ii):m_ind_vec(ii)+m_vec(ii)-1
%         eval( ['plot(tgrid,theta' num2str(i) '_opt);']);
%         if  i == m_ind_vec(ii)
%             hold on
%         end
%     end
%     xlabel('$t$','interpreter','latex');
%     ylabel('$\theta(t)$','interpreter','latex');
%     hold on
%     grid on
