clear all;
clear all;
clc;
import casadi.*
close all
%%
[settings] = default_settings_nosnoc();  
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 2;

settings.mpcc_mode = 5;
settings.opts_ipopt.ipopt.max_iter = 3e2;
settings.print_level = 2;
settings.N_homotopy = 10;
settings.initial_lambda_0 = 0; settings.initial_lambda_1 = 0; settings.initial_alpha = 0;
settings.use_fesd = 0;
settings.cross_comp_mode = 2;
settings.time_freezing = 1;

settings.pss_mode = 'Step';
settings.pss_lift_step_functions = 0;

model.mu = 0.5*1;
%%
g = 9.81;
% Symbolic variables and bounds
q = SX.sym('q',2); v = SX.sym('v',2); 
model.x = [q;v]; 
model.e = 0;
model.a_n = g;
model.x0 = [0;0.5;2;0]; 
model.f = [0;-g+g*q(1)];
model.c = q(2);
%% Simulation setings
N_finite_elements = 1;
T_sim = 1.5;
N_sim = 20;
model.T_sim = T_sim;
model.N_finite_elements = N_finite_elements;
model.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 1;
%% Call nosnoc Integrator
[results,stats,model] = integrator_fesd(model,settings);
%% read and plot results
qx = results.x_res(1,:);
qy = results.x_res(2,:);
vx = results.x_res(3,:);
vy = results.x_res(4,:);
t_opt = results.x_res(5,:);
figure
subplot(121)
plot(qx,qy);
grid on
xlabel('$q_x$','interpreter','latex');
ylabel('$q_y$','interpreter','latex');
subplot(122)
plot(t_opt,vy);
hold on
plot(t_opt,vx);
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
%%
figure
n_alpha = model.n_alpha;
for ii = 1:n_alpha
subplot(1,n_alpha,ii)
plot(results.t_grid,[results.alpha_res(ii,:),nan])
xlabel('$\tau$','interpreter','latex');
ylabel(['$\alpha_' num2str(ii) '$'],'interpreter','latex');
grid on
ylim([-0.1 1.1]);
end


%% complementarity residuals
if 0 
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
end