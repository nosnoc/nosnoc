clc; close all; clear all;
import casadi.*
%%
plot_results = false;
%% settings
% collocation settings
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
problem_options.n_s = 2;
% problem_options.time_freezing_time_rescaling = 1;
problem_options.use_speed_of_time_variables =  1; 
problem_options.local_speed_of_time_variable = 1;  
problem_options.stagewise_clock_constraint = 1;
problem_options.time_freezing = 1;
problem_options.N_stages = 40;
problem_options.N_finite_elements = 4;
% solver_options.homotopy_update_rule = 'superlinear';
solver_options.N_homotopy = 6;
solver_options.print_level = 5;
solver_options.opts_casadi_nlp.ipopt.max_iter = 1e3;
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
%% Generate Model
% angulary velocity of reference
omega = -2*pi; % no impacts
omega = -3*pi; % impacts
N_periods = 2;
model = elastic_ball_in_box_model(omega, N_periods);
problem_options.T = N_periods*(2*pi/abs(omega));
%% Solve OCP via nosnoc
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();
%% read and plot results
if plot_results
    x_res = ocp_solver.get('x');
    u_opt = ocp_solver.get('u');
    t_grid = ocp_solver.get_time_grid();
    t_grid_u = ocp_solver.get_control_grid();
    plot_results_elastic_ball_in_box
end
