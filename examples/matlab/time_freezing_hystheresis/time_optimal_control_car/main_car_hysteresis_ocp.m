%
%    This file is part of NOS-NOC.
%
%    NOS-NOC -- A software for NOnSmooth Numerical Optimal Control.
%    Copyright (C) 2022 Armin Nurkanovic, Moritz Diehl (ALU Freiburg).
%
%    NOS-NOC is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 3 of the License, or (at your option) any later version.
%
%    NOS-NOC is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with NOS-NOC; if not, write to the Free Software Foundation,
%    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%
clear all
clc
close all
import casadi.*

%% Settings
[settings] = default_settings_fesd();
settings.n_s = 2;                            % IRK stages
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
