import casadi.*
clear all
clc 
close all
%% Descrpition
% this is an optimal control example from 
% Baumrucker, Brian T., and Lorenz T. Biegler. "MPEC strategies for optimization of a class of hybrid dynamic systems." Journal of Process Control 19.8 (2009): 1248-1256.

%% settings
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
problem_options.rk_scheme = RKSchemes.GAUSS_LEGENDRE;   
problem_options.n_s = 1;                
% solver_options.N_homotopy = 8;
solver_options.homotopy_update_rule = 'superlinear';
%% Generate Model
model = tank_cascade();
%% Discretization parameters
problem_options.N_stages = 100; % number of control intervals
problem_options.N_finite_elements = 2; % number of finite element on every control interval
problem_options.T = 100;    % yime horizon
%% Solve OCP via nosnoc
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();
t_grid = ocp_solver.get_time_grid();
t_grid_u = ocp_solver.get_control_grid();
x_res = ocp_solver.get("x");
u_res = ocp_solver.get("u");
%% Read and plot Result 
figure
subplot(211)
plot(t_grid, x_res);
hold on
xlabel('$t$','interpreter','latex');
ylabel('$L(t)$','interpreter','latex');
subplot(212)
stairs(t_grid_u,[nan*ones(model.dims.n_u,1),u_res]');
hold on
xlabel('$t$','interpreter','latex');
ylabel('$u(t)$','interpreter','latex');














