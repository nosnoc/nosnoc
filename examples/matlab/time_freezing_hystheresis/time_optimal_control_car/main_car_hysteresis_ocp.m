clear all
clc
close all
import casadi.*

%% Settings
[settings] = default_settings_fesd();
settings.d = 2;                            % IRK order
settings.mpcc_mode = 6;                    % 1 - extact, 2 - smooth  ,3 -relax , 4 - penalty, 5 - elastic mode
settings.s_elastic_max = 1e0;              % upper bound for elastic variables
settings.opts_ipopt.ipopt.max_iter = 1.5e3;
settings.opts_ipopt.ipopt.print_level = 0;
settings.initial_theta = 1/2;
settings.initial_lambda = 1/2;
settings.initial_mu = 1/2;

% Step Equlibration
settings.step_equilibration  = 1;
settings.step_equilibration_mode  = 1;
settings.step_equilibration_sigma = 0.1;
settings.heuristic_step_equilibration  = 0;
settings.heuristic_step_equilibration_mode  = 1;
settings.step_equilibration_penalty = 5;
%% Time settings
settings.time_freezing = 1;
settings.time_optimal_problem = 1;
% Time freezing scaling / Speed of Time
settings.s_sot_max = 10;
settings.use_speed_of_time_variables =  1; % introduce s_tof for e
settings.local_speed_of_time_variable = 1;
% Grid settings
settings.equidistant_control_grid = 1;
settings.stagewise_clock_constraint = 1;
%% Model Settings
model.active_control = 1; % OCP or simulation problem
model.use_hystereis_model = 1;
model.linear_auxiliary_dynamics =  1; % constant or linear auxiliary dynamics
model.smooth_model = 0;
model.time_optimal_problem = settings.time_optimal_problem;
model.fuel_cost_on = 0;
%% solve OCP
[results,stats,solver_initalization,settings,model] = solve_car_hysteresis_ocp(settings,model);
%% Read and plot Result
if results.status == 1
    fprintf(results.messages.solver_message)
    fprintf(results.messages.time_message)
    fprintf(results.messages.terminal_error)
    plot_results_car_hysteresis(results,settings,model,stats)   
else
    fprintf('diverged.\n')
end
