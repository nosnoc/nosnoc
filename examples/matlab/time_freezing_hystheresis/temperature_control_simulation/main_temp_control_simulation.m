%
%    This file is part of NOSNOC.
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
% collocation settings
settings = default_settings_nosnoc();

settings.n_s = 2;                            % Degree of interpolating polynomial
% MPCC settings
settings.mpcc_mode = 3;                    % 1 - extact, 2 - smooth  ,3 -relax , 4 - penalty, 5 - elastic mode
settings.s_elastic_max = 1e1;              % upper bound for elastic variables
settings.objective_scaling_direct = 0;                   % in penalty methods  1: J = J+(1/p)*J_comp (direct)  , 0 : J = p*J+J_comp (inverse)
% Penalty/Relaxation paraemetr
settings.sigma_0 = 1e1;                     % starting smouothing parameter
settings.sigma_N = 1e-10;                   % end smoothing parameter
settings.homotopy_update_slope = 0.1;                      % decrease rate
settings.N_homotopy = ceil(abs(log(settings.sigma_N/settings.sigma_0)/log(settings.homotopy_update_slope)))+1 ;% number of steps
settings.comp_tol = 1e-14;

% time freezing settings 
settings.initial_theta = 0.5;

%^ IPOPT Settings
opts_ipopt.verbose = false;
opts_ipopt.ipopt.max_iter = 250;
opts_ipopt.ipopt.print_level = 0;
opts_ipopt.ipopt.mu_strategy = 'adaptive';
opts_ipopt.ipopt.mu_oracle = 'quality-function';
settings.opts_ipopt = opts_ipopt;
%% Generate Model
model = temp_control_model_voronoi();
%% - Simulation settings
model.T_sim = 3;
model.N_stages = 1;
model.N_finite_elements = 2;
model.N_sim = 30;
settings.use_previous_solution_as_initial_guess = 1;
%% Call FESD Integrator 
[results,stats] = integrator_fesd(model,settings);
%% Read and plot result
plot_results_for_paper
