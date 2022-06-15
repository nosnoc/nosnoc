clear all;
clear all;
clc;
import casadi.*
%%
[settings] = default_settings_nosnoc();  
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 2;

settings.mpcc_mode = 5;
settings.opts_ipopt.ipopt.max_iter = 3e2;
settings.print_level = 2;
settings.N_homotopy = 10;
settings.initial_lambda_0 = 0; settings.initial_lambda_1 = 0; settings.initial_alpha = 0;
settings.use_fesd = 1;
settings.cross_comp_mode = 2;
settings.time_freezing = 1;

settings.pss_mode = 'Step';
settings.pss_lift_step_functions= 0;

%%
% Symbolic variables and bounds
q = SX.sym('q',3); v = SX.sym('v',3); 
model.x = [q;v]; 
model.e = 0;
model.mu = 0.1;
model.a_n = 30;
model.x0 = [0;0;2;...
            3;2.5;1]; 
model.f = [0;0;-9.81];
model.c = q(3);
model.tangent1 = [1;0;0];
model.tangent2 = [0;1;0];
%% Simulation setings
N_finite_elements = 2;
T_sim = 3;
N_sim = 20;
model.T_sim = T_sim;
model.N_FE = N_finite_elements;
model.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 1;
%% Call FESD Integrator
[results,stats,model] = integrator_fesd(model,settings);
%%
qx = results.x_res(1,:);
qy = results.x_res(2,:);
qz = results.x_res(3,:);
vx = results.x_res(4,:);
vy = results.x_res(5,:);
vz = results.x_res(6,:);
t_opt = results.x_res(7,:);
figure
plot3(qx,qy,qz);
grid on
xlabel('$q_x$','Interpreter','latex');
ylabel('$q_y$','Interpreter','latex');
zlabel('$q_z$','Interpreter','latex');
%%
figure
plot(t_opt,vx);
hold on
plot(t_opt,vy);
plot(t_opt,vz);
xlabel('$t$','Interpreter','latex');
ylabel('$v$','Interpreter','latex');
% legend({'v_x','v_y','v_z'},'Interpreter','latex');
%%

%%
lambda0 = results.lambda_0_res;
lambda1 = results.lambda_1_res;
alpha = results.alpha_res;
comp1 = alpha.*lambda0;
comp2 = lambda1.*(ones(size(alpha))-alpha);
figure
subplot(121)
plot(comp1')
subplot(122)
plot(comp2')