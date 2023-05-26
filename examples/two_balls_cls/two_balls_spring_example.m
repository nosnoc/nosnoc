clear all;
clc;
import casadi.*
close all
%%
settings = NosnocOptions();
settings.irk_scheme = IRKSchemes.GAUSS_LEGENDRE;
% settings.irk_representation = 'differential';
settings.n_s = 3;
settings.print_level = 3;
% settings.N_homotopy = 8;
settings.cross_comp_mode = 1;
settings.dcs_mode = DcsMode.CLS;
settings.multiple_solvers = 0;
settings.mpcc_mode = "Scholtes_ineq";
% settings.elasticity_mode = ElasticityMode.ELL_INF;
% settings.psi_fun_type = CFunctionType.BILINEAR_TWO_SIDED;
% settings.relaxation_method = RelaxationMode.TWO_SIDED;
settings.no_initial_impacts = 1;
settings.print_details_if_infeasible = 0;
settings.pause_homotopy_solver_if_infeasible = 0;
% settings.opts_ipopt.ipopt.linear_solver = 'ma97';
settings.sigma_0 = 5;
settings.homotopy_update_slope = 0.1;
settings.real_time_plot = 0;
settings.comp_tol  = 1e-13;
settings.sigma_N = 1e-13;

settings.mpcc_mode = "elastic_ineq";
settings.elastic_scholtes = 1;
settings.sigma_0 = 1e0;
settings.homotopy_update_slope = 0.2;
%%
g = 9.81;
R = 0.2;
k = 1e4;
l = 1;
m = 1;
% Symbolic variables and bounds
q = SX.sym('q',2);
v = SX.sym('v',2);
model = NosnocModel();
model.M = eye(2);
model.x = [q;v];
model.e = 0.8;
model.mu_f = 0;
x0 = [1;2;0;0];
model.x0 = x0;
model.f_v = [-m*g+k*(q(2)-q(1)-l);-m*g-k*(q(2)-q(1)-l)];
model.f_c = q(1)-R;

%% Simulation settings
N_FE = 2;
T_sim = 1;
N_sim = 200;
model.T_sim = T_sim;
settings.N_finite_elements = N_FE;
model.N_sim = N_sim;

%% MATLAB solution
settings.use_previous_solution_as_initial_guess = 1;
[t_grid_matlab, x_traj_matlab, n_bounces, lambda_normal_guess] = two_balls_spring_matlab(T_sim,x0,model.e,1e-13);


%% Call nosnoc Integrator
initial_guess = struct();
initial_guess.x_traj = x_traj_matlab;
initial_guess.t_grid = t_grid_matlab;
initial_guess.lambda_normal_traj = lambda_normal_guess;

% [results,stats,solver] = integrator_fesd(model, settings, [], initial_guess);
[results,stats,solver] = integrator_fesd(model, settings);

%%
plot_two_ball_traj(results);

%% compare
error = norm(x_traj_matlab(end,:)'-results.x(:,end));
fprintf('Numerical error %2.2e \n',error);


