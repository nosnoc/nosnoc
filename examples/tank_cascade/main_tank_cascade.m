import casadi.*
clear all
clc 
close all
%% Descrpition
% this is an optimal control example from 
% Baumrucker, Brian T., and Lorenz T. Biegler. "MPEC strategies for optimization of a class of hybrid dynamic systems." Journal of Process Control 19.8 (2009): 1248-1256.

%% settings
settings = default_settings_nosnoc();
settings.irk_scheme = 'Gauss-Legendre';   
settings.n_s = 1;                
% settings.N_homotopy = 8;
settings.homotopy_update_rule = 'superlinear';
%% Generate Model
model = tank_cascade();
%% Discretization parameters
model.N_stages = 100; % number of control intervals
model.N_finite_elements = 2; % number of finite element on every control intevral
model.T = 100;    % yime horizon
%% Solve OCP via nosnoc
[results,stats,model,settings] = nosnoc_solver(model,settings);
%% Get variables into main workspace
unfold_struct(results,'base');
%% Read and plot Result 
figure
subplot(211)
plot(t_grid,x_opt);
hold on
xlabel('$t$','interpreter','latex');
ylabel('$L(t)$','interpreter','latex');
subplot(212)
stairs(t_grid_u,[nan*ones(n_u,1),u_opt]');
hold on
xlabel('$t$','interpreter','latex');
ylabel('$u(t)$','interpreter','latex');














