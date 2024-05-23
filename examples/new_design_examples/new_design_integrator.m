
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
problem_options.n_s = 1;
problem_options.rk_scheme = RKSchemes.GAUSS_LEGENDRE;
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.rk_representation= 'differential';
problem_options.rk_representation= 'integral';
problem_options.dcs_mode = 'Stewart';
problem_options.N_sim = N_sim;
problem_options.N_finite_elements = 2;
problem_options.T_sim = T_sim;

model = nosnoc.model.Pss();

model.x0 = -0.50;
x = SX.sym('x',1);
model.x = x;
model.c = x;
model.S = [-1; 1];
f_1 = [1]; f_2 = [-1];
model.F = [f_1 f_2];

integrator = nosnoc.Integrator(model, problem_options, solver_options);
[t_grid, x_res] = integrator.simulate();
%
figure
plot(t_grid, x_res)
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$x(t)$','Interpreter','latex')
grid on



