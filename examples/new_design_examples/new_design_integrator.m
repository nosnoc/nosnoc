
clear all
clc
close all
import casadi.*
%% discretization parameters
N_sim = 10;
T_sim = 0.75;

%% NOSNOC settings
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
problem_options.n_s = 2;
problem_options.rk_scheme = RKSchemes.GAUSS_LEGENDRE;
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.rk_representation= 'differential';
problem_options.rk_representation= 'integral';
problem_options.dcs_mode = 'Stewart';
problem_options.N_sim = N_sim;
problem_options.N_finite_elements = 2;
problem_options.T_sim = T_sim;
problem_options.use_fesd = 1;
problem_options.print_level = 5;

model = nosnoc.model.Pss();

model.x0 = -0.50;
x = SX.sym('x',1);
model.x = x;
model.c = x;
model.S = [-1; 1];
f_1 = [1]; f_2 = [-1];
model.F = [f_1 f_2];

integrator = nosnoc.integrator.FESD(model, problem_options, solver_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();
% Plotting step end points and full x vs t
figure
plot(t_grid, x_res)
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$x(t)$','Interpreter','latex')
grid on
figure
plot(t_grid_full, x_res_full)
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$x(t)$','Interpreter','latex')
grid on

% or an alternate way to get data:
t_grid = integrator.get_time_grid();
t_grid_full = integrator.get_time_grid_full();
x_res = integrator.get('x');
lam_res = integrator.get('lambda');
lam_res_full = integrator.get_full('lambda');
theta_res = integrator.get('theta');
theta_res_full = integrator.get_full('theta');
figure
plot(t_grid, x_res)
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$x(t)$','Interpreter','latex')
grid on

figure
plot(t_grid, lam_res)
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$\lambda(t)$','Interpreter','latex')
grid on
figure
plot(t_grid_full, lam_res_full)
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$\lambda(t)$','Interpreter','latex')
grid on

figure
plot(t_grid(2:end), theta_res)
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$\theta(t)$','Interpreter','latex')
grid on
figure
plot(t_grid_full(2:end), theta_res_full)
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$\theta(t)$','Interpreter','latex')
grid on
