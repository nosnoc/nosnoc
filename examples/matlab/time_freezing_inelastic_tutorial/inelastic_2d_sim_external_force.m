clear all;
clear all;
clc;
import casadi.*
close all
%%
[settings] = default_settings_nosnoc();  
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 1;

settings.mpcc_mode = 3;

settings.opts_ipopt.ipopt.max_iter = 5e2;
settings.print_level = 2;
settings.N_homotopy = 12;
settings.use_fesd = 1;
settings.cross_comp_mode = 8;
settings.time_freezing = 1;
settings.pss_lift_step_functions = 1;
settings.time_freezing_reduced_model = 0;

settings.impose_terminal_phyisical_time = 1;
settings.local_speed_of_time_variable = 1;
settings.stagewise_clock_constraint = 0;

% settings.rho_h = 0;
settings.delta_h_regularization = 1;
settings.step_equilibration = 0;
%%
g = 10;
vertical_force = 0;
% Symbolic variables and bounds
q = SX.sym('q',2); v = SX.sym('v',2); 
u = SX.sym('u',2); 
model.x = [q;v]; 
model.u = u;
model.e = 0;
model.mu = 0.3;
model.a_n = g;
model.a_n = 200;
model.x0 = [0;1;3;0]; 
model.f = [0;-g+vertical_force*g*q(1)]+u;
model.c = q(2);
model.tangent1 = [1; 0];

settings.rho_sot = 0;
%% Simulation setings
N_FE = 2;
T_sim = 1.5;
N_sim = 20;
u_sim = 1*ones(2,N_sim);
model.T_sim = T_sim;
model.N_FE = N_FE;
model.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 0;
%% Call nosnoc Integrator
[results,stats,model,settings,solver,solver_initalization] = integrator_fesd(model,settings,u_sim);
if 0
    % re-run simulation without creating new solver object
    [results_rerun,stats,model] = integrator_fesd(model,settings,u_sim,solver,solver_initalization);
    norm(results.x_res-results_rerun.x_res)
end
%% read and plot results
unfold_struct(results_rerun,'base');

qx = x_res(1,:);
qy = x_res(2,:);
vx = x_res(3,:);
vy = x_res(4,:);
t_opt = x_res(5,:);
figure
subplot(121)
plot(qx,qy);
axis equal
grid on
xlabel('$q_x$','interpreter','latex');
ylabel('$q_y$','interpreter','latex');
% axis equal
subplot(122)
plot(t_opt,vy);
hold on
plot(t_opt,vx);
% axis equal
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');

% [qx(end),qy(end),t_opt(end)]

%%
alpha1 = alpha_res(1,:);
alpha2 = alpha_res(2,:);
alpha3 = alpha_res(3,:);
theta1 = alpha1+(1-alpha1).*(alpha2);
alpha_aux = (1-alpha1).*(1-alpha2);
theta2 = alpha_aux.*(1-alpha3);
theta3 = alpha_aux.*(alpha3);
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
%% speed of time
figure
subplot(121)
plot(t_grid,t_opt)
hold on
plot(t_grid,t_grid,'k--')
grid on
xlabel('$\tau$','interpreter','latex');
ylabel('$t$','interpreter','latex');
subplot(122)
stairs(s_sot_res)
grid on
xlabel('simulation step','interpreter','latex');
ylabel('$s$','interpreter','latex');



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