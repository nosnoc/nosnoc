clear all;
clear all;
clc;
import casadi.*
close all
%%
[settings] = default_settings_nosnoc();  
settings.irk_scheme = 'Radau-IIA';
% settings.irk_scheme = 'Gauss-Legendre';
settings.n_s = 3;
settings.mpcc_mode = 2;
settings.print_level = 2;
settings.use_fesd = 1;
settings.comp_tol = 1e-7;
settings.time_freezing = 1;
settings.pss_lift_step_functions = 1;

settings.impose_terminal_phyisical_time = 1;
settings.local_speed_of_time_variable = 1;
settings.stagewise_clock_constraint = 0;

% settings.rho_h = 0;
% settings.delta_h_regularization = 1;
% settings.step_equilibration = 0;
%%
g = 10;
% Symbolic variables and bounds
q = SX.sym('q',1); 
v = SX.sym('v',1); 
model.x = [q;v]; 
model.e = 0.9;
model.k_aux = 50;
model.x0 = [1;0]; 
model.f = -g;
model.c = q;
%% Simulation setings
N_FE = 3;
T_sim = 5;
N_sim = 50;
model.T_sim = T_sim;
model.N_FE = N_FE;
model.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 0;
%% Call nosnoc Integrator
[results,stats,model] = integrator_fesd(model,settings);
%% read and plot results
unfold_struct(results,'base');
qx = x_res(1,:);
vx = x_res(2,:);
t_opt = x_res(3,:);
figure
subplot(121)
plot(t_opt,qx);
axis equal
grid on
xlabel('$t$','interpreter','latex');
ylabel('$q$','interpreter','latex');
% axis equal
subplot(122)
plot(t_opt,vx);
hold on
% axis equal
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');

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
