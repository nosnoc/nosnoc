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

model.mu = 0.2*1;
%%
% Symbolic variables and bounds
q = SX.sym('q',2); v = SX.sym('v',2); 
model.x = [q;v]; 
model.e = 0;

model.a_n = 10;
model.x0 = [5;0.5;1;1]; 
model.f = [0;-9.81];
model.c = q(2);
%% Simulation setings
N_finite_elements = 2;
T_sim = 1.5;
N_sim = 20;
model.T_sim = T_sim;
model.N_finite_elements = N_finite_elements;
model.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 1;
%% Call FESD Integrator
[results,stats,model] = integrator_fesd(model,settings);
%%
qx = results.x_res(1,:);
qy = results.x_res(2,:);
vx = results.x_res(3,:);
vy = results.x_res(4,:);
t_opt = results.x_res(5,:);
figure
subplot(121)
plot(qx,qy);
subplot(122)
plot(t_opt,vy);
hold on
plot(t_opt,vx);

%% complementarity residuals
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