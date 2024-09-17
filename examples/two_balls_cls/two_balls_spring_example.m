clear all;
clc;
import casadi.*
close all
%%
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
% problem_options.rk_representation = 'differential';
problem_options.n_s = 3;
solver_options.print_level = 3;
% solver_options.N_homotopy = 8;
problem_options.cross_comp_mode = 'FE_FE';
problem_options.dcs_mode = DcsMode.CLS;
problem_options.no_initial_impacts = 1;
%problem_options.relax_terminal_numerical_time = 'ELL_1';
%problem_options.rho_terminal_numerical_time = 1e5;

solver_options.homotopy_steering_strategy = HomotopySteeringStrategy.ELL_INF;
solver_options.print_details_if_infeasible = 0;
solver_options.pause_homotopy_solver_if_infeasible = 0;
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
solver_options.real_time_plot = 0;
solver_options.complementarity_tol = 1e-10;
solver_options.sigma_N = 1e-11;
solver_options.decreasing_s_elastic_upper_bound = 1;
solver_options.sigma_0 = 5;
solver_options.homotopy_update_slope = 0.1;
solver_options.opts_casadi_nlp.ipopt.max_iter = 5e3;

%%
g = 9.81;
R = 0.2;
k = 1e4;
l = 1;
m = 1;
% Symbolic variables and bounds
q = SX.sym('q',2);
v = SX.sym('v',2);
model = nosnoc.model.Cls();
model.M = eye(2);
model.x = [q;v];
model.e = 0.8;
model.mu = 0;
x0 = [1;2;0;0];
model.x0 = x0;
model.f_v = [-m*g+k*(q(2)-q(1)-l);-m*g-k*(q(2)-q(1)-l)];
model.f_c = q(1)-R;

%% Simulation settings
N_FE = 2;
T_sim = 1;
N_sim = 200;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;

%% MATLAB solution
solver_options.use_previous_solution_as_initial_guess = 0;
[t_grid_matlab, x_traj_matlab, n_bounces, lambda_normal_guess] = two_balls_spring_matlab(T_sim,x0,model.e,1e-13);


%% Call nosnoc Integrator
initial_guess = struct();
initial_guess.x_traj = x_traj_matlab;
initial_guess.t_grid = t_grid_matlab;
initial_guess.lambda_normal_traj = lambda_normal_guess;

% [results,stats,solver] = integrator_fesd(model, settings, [], initial_guess);
integrator = nosnoc.Integrator(model, problem_options, solver_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();
%%
plot_two_ball_traj(results);

%% compare
error = norm(x_traj_matlab(end,:)'-results.x(:,end));
fprintf('Numerical error %2.2e \n',error);


