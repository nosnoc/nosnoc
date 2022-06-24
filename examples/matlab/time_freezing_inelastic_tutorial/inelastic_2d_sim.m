clear all;
clear all;
clc;
import casadi.*
close all
%%
[settings] = default_settings_nosnoc();  
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 2;

settings.mpcc_mode = 3;

settings.opts_ipopt.ipopt.max_iter = 5e2;
settings.print_level = 2;
settings.N_homotopy = 12;
settings.use_fesd = 1;
settings.cross_comp_mode = 8;
settings.time_freezing = 1;
settings.pss_lift_step_functions = 1;
settings.time_freezing_reduced_model = 0;

model.mu = 0.3*1;
%%
g = 9.81;
vertical_force = 0;
% Symbolic variables and bounds
q = SX.sym('q',2); v = SX.sym('v',2); 
model.x = [q;v]; 
model.e = 0;
model.a_n = g;
model.x0 = [0;0.2;1;0]; 
model.x0 = [0;0.0;1;0]; 
model.f = [0;-g+vertical_force*g*q(1)];
model.c = q(2);
%% Simulation setings
N_FE = 2;
T_sim = 1.5;
N_sim = 20;
model.T_sim = T_sim;
model.N_FE = N_FE;
model.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 0;
%% Call nosnoc Integrator
[results,stats,model] = integrator_fesd(model,settings);
%% read and plot results
unfold_struct(results,'base');
qx = x_res(1,:);
qy = x_res(2,:);
vx = x_res(3,:);
vy = x_res(4,:);
t_opt = x_res(5,:);
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

[qx(end),qy(end),t_opt(end)]
%%
figure
n_alpha = model.n_alpha;
for ii = 1:n_alpha
subplot(1,n_alpha,ii)
plot(t_grid,[alpha_res(ii,:),nan])
xlabel('$\tau$','interpreter','latex');
ylabel(['$\alpha_' num2str(ii) '$'],'interpreter','latex');
grid on
ylim([-0.1 1.1]);
end
%%
alpha1 = alpha_res(1,:);
alpha2 = alpha_res(2,:);
alpha3 = alpha_res(3,:);
theta1 = alpha1+(1-alpha1).*(alpha2);
theta2 = (1-alpha1).*(1-alpha2).*(1-alpha3);
theta3 = (1-alpha1).*(1-alpha2).*(alpha3);
figure;
subplot(131)
plot(t_grid,[theta1,nan])
xlabel('$\tau$','interpreter','latex');
ylabel(['$\theta_1$'],'interpreter','latex');
grid on
ylim([-0.1 1.1]);
subplot(132)
plot(t_grid,[theta2,nan])
xlabel('$\tau$','interpreter','latex');
ylabel(['$\theta_2$'],'interpreter','latex');
grid on
ylim([-0.1 1.1]);
subplot(133)
plot(t_grid,[theta3,nan])
xlabel('$\tau$','interpreter','latex');
ylabel(['$\theta_3$'],'interpreter','latex');
grid on
ylim([-0.1 1.1]);


%% complementarity residuals
if 0 
lambda0 = lambda_0_res;
lambda1 = lambda_1_res;
alpha = alpha_res;
comp1 = alpha.*lambda0;
comp2 = lambda1.*(ones(size(alpha))-alpha);
figure
subplot(121)
plot(comp1')
subplot(122)
plot(comp2')
end