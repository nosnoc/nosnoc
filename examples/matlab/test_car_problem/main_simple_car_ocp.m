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

%% settings
[settings] = default_settings_fesd();  %% Optionally call this function to have an overview of all options. 
                                        % Missing settings are anyway filled in latter with their respecitve values.
% irk settings
settings.n_s = 2;                            % Degree of interpolating polynomial
settings.irk_scheme= 'radau';     % Collocation scheme: radau or legendre

% MPCC settings
settings.mpcc_mode = 5;                    % 1 - extact, 2 - smooth  ,3 -relax , 4 - penalty, 5 - elastic mode
settings.s_elastic_max = 1e2;              % upper bound for elastic variables
settings.s_elastic_min = 0;              % upper bound for elastic variables
settings.s_elastic_0 = 1;              % upper bound for elastic variables
settings.objective_scaling_direct = 1;      % in penalty methods  1: J = J+(1/p)*J_comp (direct)  , 0 : J = p*J+J_comp (inverse)
settings.lp_initalization = 0;
settings.initial_theta = 0.5*0;
settings.initial_lambda = 0.5*0;
settings.initial_mu = 0.5*0;
% Penalty/Relaxation paraemetr
settings.sigma_0 = 1e2;                     % starting smouothing parameter
settings.sigma_N = 1e-14;                   % end smoothing parameter
settings.kappa = 0.1;                      % decrease rate
settings.N_homotopy = ceil(abs(log(settings.sigma_N/settings.sigma_0)/log(settings.kappa)))+1 ;% number of steps
settings.N_homotopy = 8;

% settings.N_homotopy = 1;% number of steps
settings.comp_tol = 1e-14;
settings.cross_comp_mode = 3;

%^ IPOPT Settings
opts_ipopt.verbose = false;
opts_ipopt.ipopt.max_iter = 600;
opts_ipopt.ipopt.max_iter = 1000;
opts_ipopt.ipopt.print_level = 5;
opts_ipopt.ipopt.mu_strategy = 'adaptive';
opts_ipopt.ipopt.mu_oracle = 'quality-function';
settings.opts_ipopt = opts_ipopt;

% finite elements with switch detection
settings.use_fesd = 1;       % turn on moving finite elements algortihm
settings.store_all_homotopy_iterates = 1;
% step equilibration
settings.gamma_h = 1;                    % how much can h addapt
% regularize_h --> heuristic_step_equilibration 
settings.regularize_h = 1;                 % add quadratic cost term rho_h = (h-hk)^2
settings.rho_h = 1;                        % regularization penalty
settings.heuristic_step_equilibration = 0;
settings.step_equilibration = 1;
settings.step_equilibration_mode = 1;
settings.step_equilibration_sigma = 0.01;
settings.treat_step_equilibration_with_mpcc_method = 0;
settings.step_equilibration_penalty = 1e2;

settings.convex_sigma_rho_constraint = 1;
settings.sigma_penalty = 1e2;
settings.rho_penalty = 1;

%% Time settings

% Tf % Optimal control
settings.time_freezing = 0;
settings.time_optimal_problem = 1;

% Time freezing scaling 
settings.use_speed_of_time_variables =  1; % introduce s_tof for e
settings.local_speed_of_time_variable = 1;  

% Grid settings
settings.equidistant_control_grid = 1;
settings.stagewise_clock_constraint = 0;
[results,stats,solver_initalization,settings,model] = solve_simple_car_ocp(settings);

%% Read and plot Result
if results.status == 1
fprintf(results.messages.solver_message)
fprintf(results.messages.time_message)
fprintf(results.messages.objective_message)
plot_results_simple_car(results,settings,model,stats)
else
    fprintf('diverged.\n')
end
% %%
% figure
% stairs(results.w_opt(model.ind_h))
% grid on

