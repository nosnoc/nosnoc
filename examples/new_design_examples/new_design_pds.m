
clear all
clc
close all
import casadi.*
%% discretization parameters
N_sim = 100;
T_sim = 10;

%% NOSNOC settings
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
problem_options.n_s = 2;
settings.homotopy_update_slope = 0.1;
%problem_options.rk_scheme = RKSchemes.GAUSS_LEGENDRE;
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.rk_representation= 'differential_lift_x'; 
%problem_options.rk_representation= 'integral';
problem_options.N_sim = N_sim;
problem_options.N_finite_elements = 2;
problem_options.T_sim = T_sim;
solver_options.homotopy_steering_strategy = 'ELL_INF';


model = nosnoc.model.Pds();

model.x0 = [sqrt(2)/2;sqrt(2)/2];
x = SX.sym('x',2);
model.x = x;
model.c = x(2)+0.2;
model.f_x = [x(2);-x(1)];

model.x0 = [0;pi-2;0];

integrator = nosnoc.Integrator(model, problem_options, solver_options);
[t_grid, x_res] = integrator.simulate();
%
figure
plot(x_res(1,:), x_res(2,:))
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$x(t)$','Interpreter','latex')
grid on

