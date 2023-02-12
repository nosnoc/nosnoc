clear all;
close all;
clc;
import casadi.*
%%
[settings] = default_settings_nosnoc();  
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 2;
settings.N_homotopy = 6;
settings.homotopy_update_rule = 'superlinear';
settings.print_level = 3;
settings.time_freezing = 1;
settings.s_sot_max = 10;
settings.s_sot_min = 0.1;
settings.equidistant_control_grid = 1;
settings.pss_lift_step_functions = 1;
settings.polishing_step = 1;
% settings.opts_ipopt.ipopt.linear_solver = 'ma57';
%%
g = 9.81;
u_max = 10;
% Symbolic variables and bounds
q = SX.sym('q',2);
v = SX.sym('v',2); 
x = [q;v];
u = SX.sym('u');
N_FE = 3;
model.T = 4;
model.N_stages = 15;
model.N_finite_elements  = N_FE;
model.x = x;
model.u = u;
model.e = 0;
model.mu = 0.7;
model.a_n = g;
model.x0 = [0;1;0;0]; 
model.f_v = [0+u;-g];
model.J_tangent = [1;0];
model.c = q(2);
model.g_terminal = [x-[2;0;0;0]];
model.lbu = -u_max;
model.ubu = u_max;
model.f_q = u'*u;
%% Call nosnoc solver
[results,stats,model] = nosnoc_solver(model,settings);
%%
unfold_struct(results,'base')
qx = x_opt(1,:);
qy = x_opt(2,:);
vx = x_opt(3,:);
vy = x_opt(4,:);
t_opt = x_opt(5,:);
u_opt = u_opt;
s_opt = w_opt(model.ind_sot);

figure
subplot(131)
plot(qx,qy,'LineWidth',2);
xlim([0 2]);
ylim([-0.1 1]);
xlabel('$q_1$','Interpreter','latex');
ylabel('$q_2$','Interpreter','latex');
grid on
subplot(132)
plot(t_opt,vy,'LineWidth',2);
hold on
plot(t_opt,vx,'LineWidth',2);
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$v$','Interpreter','latex');
legend({'$v_1$','$v_2$'},'Interpreter','latex','Location','best');
subplot(133)
stairs(t_opt(1:N_FE:end),[u_opt,[nan]]','LineWidth',2);
hold on
stairs(t_opt(1:N_FE:end),[s_opt',nan],'LineWidth',2);
legend({'$u_1(t)$','$s(t)$'},'Interpreter','latex','location','best');
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');

%%
if 0

figure
subplot(131)
plot(alpha_opt(1,:))
grid on
ylabel('$\alpha_1$','Interpreter','latex');
subplot(132)
plot(alpha_opt(2,:))
ylabel('$\alpha_2$','Interpreter','latex');
grid on
subplot(133)
plot(alpha_opt(3,:))
ylabel('$\alpha_3$','Interpreter','latex');
grid on

theta1 = alpha_opt(1,:) + (1-alpha_opt(1,:)).*(alpha_opt(2,:));
theta2 = (1-alpha_opt(1,:)).*(1-alpha_opt(2,:)).*(1-alpha_opt(3,:));
theta3 = (1-alpha_opt(1,:)).*(1-alpha_opt(2,:)).*(alpha_opt(3,:));

figure
subplot(131)
plot(theta1)
grid on
ylabel('$\theta_1$','Interpreter','latex');
subplot(132)
plot(theta2)
ylabel('$\theta_2$','Interpreter','latex');
grid on
subplot(133)
plot(theta3)
ylabel('$\theta_3$','Interpreter','latex');
grid on

figure
subplot(131)
plot(t_grid,qy)
subplot(132)
grid on
ylabel('$f_c(q)$','Interpreter','latex');
plot(t_grid,vy)
ylabel('$n^top v$','Interpreter','latex');
grid on
subplot(133)
plot(t_grid,vx)
ylabel('$t^top v$','Interpreter','latex');
grid on

%
figure
subplot(131)
plot(lambda_1_opt(1,:))
grid on
ylabel('$\lambda^+_1$','Interpreter','latex');
subplot(132)
plot(lambda_1_opt(2,:))
ylabel('$\lambda^+2$','Interpreter','latex');
grid on
subplot(133)
plot(lambda_1_opt(3,:))
ylabel('$\lambda^+_3$','Interpreter','latex');
grid on

figure
subplot(131)
plot(lambda_0_opt(1,:))
grid on
ylabel('$\lambda^+1$','Interpreter','latex');
subplot(132)
plot(lambda_0_opt(2,:))
ylabel('$\lambda^+_2$','Interpreter','latex');
grid on
subplot(133)
plot(lambda_0_opt(3,:))
ylabel('$\lambda^+3$','Interpreter','latex');
grid on
end
