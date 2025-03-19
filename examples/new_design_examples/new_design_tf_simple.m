%clear all;
clc;
import casadi.*
close all
%% init nosnoc
problem_options = nosnoc.Options();
solver_options = nosnoc.reg_homotopy.Options();
model = nosnoc.model.Cls();
%% Simulation settings
problem_options.rk_scheme = RKSchemes.GAUSS_LEGENDRE;
%problem_options.rk_scheme = RKSchemes.RADAU_IIA;
N_FE = 3;
T_sim = 1.1;
N_sim = 21;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;
problem_options.n_s = 2;
problem_options.cross_comp_mode = 1;

% Time-freezing settings
problem_options.time_freezing = 1;
problem_options.impose_terminal_phyisical_time = 1;
problem_options.stagewise_clock_constraint = 0;
problem_options.s_sot_max = 100;
%problem_options.a_n = 100;
problem_options.k_aux = 50;
%problem_options.stabilizing_q_dynamics = true;


%problem_options.rho_terminal_physical_time = 1e4;
%problem_options.relax_terminal_physical_time = "ELL_1";

% problem_options.relax_terminal_numerical_time = ConstraintRelaxationMode.ELL_1;
% problem_options.rho_terminal_numerical_time = 1e3;

problem_options.dcs_mode = DcsMode.Heaviside;
problem_options.time_freezing_Heaviside_lifting = true;
%problem_options.rho_h = 0;
%% mpcc solver options
solver_options.print_level = 3;
solver_options.complementarity_tol = 1e-10;
solver_options.sigma_N = 1e-10;
solver_options.N_homotopy = 16;
solver_options.homotopy_steering_strategy = "ELL_INF";
solver_options.decreasing_s_elastic_upper_bound = true;
solver_options.pause_homotopy_solver_if_infeasible = 0;
solver_opts.use_previous_solution_as_initial_guess = 1;
% solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
%% Model
g = 10;
vertical_force = 0;
q = SX.sym('q',2); 
v = SX.sym('v',2);
model.x = [q;v]; 
model.e = 1*0;
model.mu = 0.0;
model.x0 = [0;1;1;0]; 
model.f_v = [0;-g+vertical_force*g*q(1)];
model.f_c = q(2);
model.J_tangent = [1; 0];
model.dims.n_dim_contact = 2;

integrator = nosnoc.Integrator(model, problem_options, solver_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();

%% plot algebraics
figure
if problem_options.dcs_mode == DcsMode.Stewart
    theta = integrator.get_full('theta');
    subplot(121)
    plot(t_grid_full, theta)
    subplot(122)
    theta = integrator.get('theta');
    plot(t_grid, theta)
else
    alpha= integrator.get_full('alpha');
    subplot(211)
    t_grid_alpha = t_grid_full;
    t_grid_alpha((problem_options.n_s+2):(problem_options.n_s+1):end) = [];
    plot(t_grid_alpha, alpha)
    subplot(212)
    alpha= integrator.get('alpha');
    plot(t_grid, alpha)

    figure
    z= integrator.get_full('z');
    subplot(211)
    plot(t_grid_alpha, z)
    subplot(212)
    z= integrator.get('z');
    plot(t_grid, z)
end

%% Plot
figure
plot(x_res(end,:)', x_res(1:2,:)')
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$x(t)$','Interpreter','latex')
figure
plot(x_res_full(end,:)', x_res_full(1:2,:)')
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$x(t)$','Interpreter','latex')
figure
plot(x_res(end,:), x_res(3:4,:))
xlabel('$t$','Interpreter','latex')
ylabel('$v(t)$','Interpreter','latex')
grid on
figure
plot(x_res_full(end,:), x_res_full(3:4,:))
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$v(t)$','Interpreter','latex')
title("full")
grid on
figure
plot(x_res_full(1,:), x_res_full(2,:))
figure
plot(x_res_full(2,:), x_res_full(4,:))
figure
hold on
plot(t_grid, x_res(end,:))
plot(t_grid, t_grid, 'k--')
hold off
