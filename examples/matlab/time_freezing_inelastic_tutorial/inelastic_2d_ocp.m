clear all;
close all;
clc;
import casadi.*
%%
[settings] = default_settings_nosnoc();  
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 2;
settings.sigma_0 = 1;
settings.mpcc_mode = 3;
settings.kappa = 0.1;
settings.N_homotopy = 6;
settings.cross_comp_mode = 3;
settings.opts_ipopt.ipopt.max_iter = 1e3;
settings.print_level = 3;

% settings.opts_ipopt.ipopt.tol = 1e-14;

settings.time_freezing = 1;
settings.s_sot_max = 10;
settings.s_sot_min = 0.1;
settings.equidistant_control_grid = 1;
settings.pss_lift_step_functions = 1;


settings.opts_ipopt.ipopt.linear_solver = 'ma57';
settings.step_equilibration = 1;
settings.step_equilibration_mode = 3;

settings.polishing_step = 1;
%%
g = 9.81;
u_max = 10;
N_stg = 20;  N_FE = 2; 
% Symbolic variables and bounds
q = SX.sym('q',2); v = SX.sym('v',2); 
x = [q;v];
u = SX.sym('u');
model.T = 2;
model.N_stages = N_stg;
model.N_finite_elements  = N_FE;
model.x = x;
model.u = u;
model.e = 0;
model.mu = 0.6;
model.a_n = g;
model.x0 = [0;1.5;0;0]; 
model.f = [0+u;-g];

model.c = q(2);
model.g_terminal = [x-[3;0;0;0]];
model.lbu = -u_max;
model.ubu = u_max;
model.f_q = u'*u;
%% Call nosnoc solver
[results,stats,model] = nosnoc_solver(model,settings);
%%
qx = results.x_opt(1,:);
qy = results.x_opt(2,:);
vx = results.x_opt(3,:);
vy = results.x_opt(4,:);
t_opt = results.x_opt(5,:);
u_opt = results.u_opt;
s_opt = results.w_opt(model.ind_sot);

figure
subplot(131)
plot(qx,qy,'LineWidth',2);
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
legend({'$v_1$','$v_2$'},'Interpreter','latex');
subplot(133)
stairs(t_opt(1:N_FE:end),[u_opt,[nan]]','LineWidth',2);
hold on
stairs(t_opt(1:N_FE:end),[s_opt',nan],'LineWidth',2);
legend({'$u_1(t)$','$s(t)$'},'Interpreter','latex','location','best');
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
%%
% saveas(gcf,'ocp_example','epsc')
