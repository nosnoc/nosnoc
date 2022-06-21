clear all;
close all;
clc;
import casadi.*
%%
[settings] = default_settings_nosnoc();  
settings.irk_scheme = 'legendre';
settings.irk_scheme = 'radau';
settings.n_s = 1;
settings.pss_mode = 'Step';
% settings.sigma_0 = 100;
settings.mpcc_mode = 2;
settings.N_homotopy = 10;
settings.cross_comp_mode = 3;
settings.opts_ipopt.ipopt.max_iter = 1e3;
settings.print_level = 3;


settings.use_fesd = 1;
settings.time_freezing = 1;
settings.stagewise_clock_constraint = 1;

settings.s_sot_max = 3;
settings.s_sot_min = 0.5;
settings.equidistant_control_grid = 1;

settings.pss_lift_step_functions = 0;


settings.opts_ipopt.ipopt.linear_solver = 'ma27';
settings.opts_ipopt.ipopt.linear_solver = 'ma57';
%%
g = 9.81;
N_stages = 15; N_finite_elements  = 3; 
% Symbolic variables and bounds
q = SX.sym('q',2); v = SX.sym('v',2); 
x = [q;v];
u = SX.sym('u',2);
model.T = 4;
model.N_stages = N_stages;
model.N_finite_elements  = N_finite_elements;
model.x = x;
model.M = eye(2);
model.u = u;
model.e = 0;
model.mu = 0.25;
model.a_n = 9.81;
model.x0 = [0;1;1;-2]; 
model.f = [0;-g]+u;
model.c = q(2);
% model.g_terminal = [x-[3;0;0;0]];
model.g_terminal = [q-[3;1]];
settings.terminal_constraint_relxataion = 0;
settings.rho_terminal = 1e3;
model.lbu = -2*g*ones(2,1);
model.ubu = 2*g*ones(2,1);
model.f_q = u'*u;
model.f_q_T = 100*v'*v;
%% Call FESD Integrator
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
stairs(t_opt(1:N_finite_elements:end),[u_opt,[nan;nan]]','LineWidth',2);
hold on
stairs(t_opt(1:N_finite_elements:end),[s_opt',nan],'LineWidth',2);
legend({'$u_1(t)$','$s(t)$'},'Interpreter','latex','location','best');
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
%%
% saveas(gcf,'ocp_example','epsc')
