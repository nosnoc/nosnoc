clear all;
clear all;
clc;
import casadi.*
%%
[settings] = default_settings_nosnoc();  
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 2;
settings.pss_mode = 'Step';
settings.pss_lift_step_functions= 1;
settings.mpcc_mode = 3;
settings.opts_ipopt.ipopt.max_iter = 3e2;
settings.print_level = 3;
settings.N_homotopy = 10;
settings.initial_lambda_0 = 0; settings.initial_lambda_1 = 0; settings.initial_alpha = 0;
settings.use_fesd = 1;
settings.cross_comp_mode = 10;
settings.time_freezing = 1;
settings.stagewise_clock_constraint = 1;


%%
N_stages = 20; N_finite_elements  = 3; 
u_max = 10;
% Symbolic variables and bounds
q = SX.sym('q',2); v = SX.sym('v',2); 
u = SX.sym('u');
model.T = 1.5;
model.N_stages = N_stages;
model.N_finite_elements  = N_finite_elements;
model.x = [q;v];
model.M = eye(2);
model.u = u;
model.e = 0;
model.mu = 0.9;
model.a_n = 30;
model.x0 = [0;0.5;1;1]; 
model.f = [u;-9.81];
model.c = q(2);
model.lbu = -u_max;
model.ubu = u_max;
model.g_terminal = [q-[3;0]];
model.f_q = u^2;
%% Call FESD Integrator
[results,stats,model] = nosnoc_solver(model,settings);
%%
qx = results.x_opt(1,:);
qy = results.x_opt(2,:);
vx = results.x_opt(3,:);
vy = results.x_opt(4,:);
t_opt = results.x_opt(5,:);
u_opt = results.u_opt;
figure
subplot(131)
plot(qx,qy);
xlabel('$q_x$','Interpreter','latex');
ylabel('$q_y$','Interpreter','latex');
grid on
subplot(132)
plot(t_opt,vy);
hold on
plot(t_opt,vx);
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$v$','Interpreter','latex');
legend({'$v_x$','$v_y$'},'Interpreter','latex');
subplot(133)
stairs(t_opt(1:N_finite_elements:end),[u_opt,nan]);
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$v$','Interpreter','latex');
%%
