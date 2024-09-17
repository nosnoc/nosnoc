clear all;
clc;
import casadi.*
close all
%% init
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
model = nosnoc.model.Cls();
%% settings
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 2;
problem_options.cross_comp_mode = 1;
problem_options.dcs_mode = DcsMode.CLS;
problem_options.no_initial_impacts = 1;
problem_options.friction_model = "Polyhedral";
%problem_options.friction_model = "Conic"; % "Conic"
%problem_options.conic_model_switch_handling = "Lp";

% Initialization
% problem_options.initial_Lambda_normal = 10;
% problem_options.initial_lambda_normal = 0;
% problem_options.initial_Y_gap = 0;
% problem_options.initial_y_gap = 0;

problem_options.fixed_eps_cls = 1;
%problem_options.relax_terminal_numerical_time = ConstraintRelaxationMode.ELL_2;
%problem_options.rho_terminal_numerical_time = 1e9;
%problem_options.relax_fesdj_impulse = ConstraintRelaxationMode.ELL_2;
%problem_options.rho_fesdj_impulse = 1e9;
problem_options.gamma_h = 1;

solver_options.complementarity_tol = 1e-6;
solver_options.use_previous_solution_as_initial_guess = 0;
solver_options.print_level = 3;
solver_options.N_homotopy = 10;
solver_options.sigma_0 = 10;
solver_options.homotopy_steering_strategy = 'ELL_INF';
solver_options.decreasing_s_elastic_upper_bound = true;
%%
g = 9.81;
% Symbolic variables and bounds
q = SX.sym('q',2);
v = SX.sym('v',2);
model.M = diag([1,1]);
model.x = [q;v];
model.e = 0;
model.mu = 0.1;
% impact
model.x0 = [0;1;4;0];
% Free flight
% model.x0 = [0;5.5;4;0];
% Friction slip slip
% model.x0 = [0;0.0;-0.6;0];
% friction slip stick
% model.x0 = [0;0.0;+0.6;0];
model.f_v = [3*0;-g];
model.f_c = q(2);
model.J_tangent = [1; 0];
model.D_tangent = [1,-1;0,0];

%% Simulation settings
N_FE = 2;
T_sim = 1.7;
%T_sim = 1.7/10;
N_sim = 10;
%N_sim = 1;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;

%% Call nosnoc Integrator
integrator = nosnoc.Integrator(model, problem_options, solver_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();
%% read and plot results
qx = x_res(1,:);
qy = x_res(2,:);
vx = x_res(3,:);
vy = x_res(4,:);
figure(1)
subplot(121)
plot(qx,qy);
axis equal
grid on
xlabel('$q_x$','interpreter','latex');
ylabel('$q_y$','interpreter','latex');
% axis equal
subplot(122)
plot(t_grid,vy);
hold on
plot(t_grid,vx);
% axis equal
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
%% forces
if N_sim == 1
figure
subplot(311)
plot(t_grid,x_res(3:4,:));
hold on
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
subplot(312)
plot(t_grid,[nan*ones(model.n_contacts,1),results.lambda_normal])
hold on
plot(t_grid,[nan*ones(model.n_tangents,1),results.lambda_tangent])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\lambda$','interpreter','latex');
subplot(313)
stem(t_grid,[results.Lambda_normal,nan*ones(model.n_contacts,1)])
hold on
stem(t_grid,[results.Lambda_tangent,nan*ones(model.n_tangents,1)]')
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda$','interpreter','latex');
xlim([0 T_sim])
end
