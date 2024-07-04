clear all
clc
close all
import casadi.*
%% discretization parameters
N_sim = 31;
T_sim = 10;

%% NOSNOC settings
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
problem_options.n_s = 3;
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
%problem_options.rk_representation= RKRepresentation.differential_lift_x; 
problem_options.rk_representation = RKRepresentation.integral;
problem_options.cross_comp_mode = CrossCompMode.STAGE_STAGE;
%problem_options.cross_comp_mode = CrossCompMode.FE_FE;
problem_options.N_sim = N_sim;
problem_options.N_finite_elements = 2;
problem_options.T_sim = T_sim;
problem_options.rho_h = 1e-4;
problem_options.gcs_lift_gap_functions = true;
problem_options.step_equilibration = StepEquilibrationMode.l2_relaxed_scaled;
%solver_options.homotopy_steering_strategy = 'ELL_INF';
solver_options.complementarity_tol = 1e-10;
solver_options.print_level = 3;

model = nosnoc.model.Pds();

model.x0 = [sqrt(2)/2;sqrt(2)/2];
x = SX.sym('x',2);
model.x = x;
model.c = x(2)+0.2;
model.f_x_unconstrained = [x(2);-x(1)];

model.x0 = [0;pi-2];

integrator = nosnoc.Integrator(model, problem_options, solver_options);
[t_grid, x_res] = integrator.simulate();
%
figure
plot(x_res(1,:), x_res(2,:))
grid on
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
grid on

c_fun = casadi.Function('c', {model.x}, {model.c});
c = full(c_fun(integrator.get_full('x')))';
x_full = integrator.get_full('x');
lambda = integrator.get('lambda');

% extened values over all rk stage points
t_grid_full = integrator.get_time_grid_full();
lambda_full = integrator.get_full('lambda');

%

figure
plot(t_grid,lambda,'LineWidth',2)
hold on
plot(t_grid,x_res,'LineWidth',1.5)
grid on
xlabel('$t$','Interpreter','latex')
legend({'$\lambda(t)$','$x_1(t)$','$x_2(t)$'},'interpreter','latex');

