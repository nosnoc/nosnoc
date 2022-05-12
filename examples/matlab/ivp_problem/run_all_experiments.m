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
%%
clear all
clc
close all
import casadi.*
%% Experiments for derivatives of objective and local minima
run_objective_function_experiment = 0;
run_local_minima_experiment = 1;


N_samples = 150;
scenario.N_samples = N_samples;
scenario.x0_vec = linspace(-2,-0.8,N_samples);
scenario.save_results = 1;
%% settings
settings = default_settings_fesd();
settings.n_s = 2;
settings.irk_scheme = 'Gauss-Legendre';
settings.equidistant_control_grid = 0;
%% Generate Model
model.T = 2;
model.N_stages = 1;
model.N_finite_elements = 25;
model.x0 = -1;
% Variable defintion
x = MX.sym('x');
model.x = x;
model.c = x;
model.S = [1;-1];
% modes of the ODE
f_11 = [1];
f_12 = [3];
model.F = [f_11 f_12];
% objective
model.f_q = x^2;
model.f_q_T = (x-5/3)^2;
if run_local_minima_experiment
    % Scenario-1  : use FESD with fixed parameters
%     scenario.scenario_name = 'local_min_GL4_fesd_fixed';
%     settings.use_fesd = 1;
%     settings.mpcc_mode = 3;
%     settings.sigma_0 = 1e-15;
%     settings.sigma_N = 1e-15;
%     settings.N_homotopy = 1;
%     [results] = local_minima_experiment(scenario,settings,model);
%     % Scenario-2  : use Standard with fixed parameters
%     scenario.scenario_name = 'local_min_GL4_std_fixed';
%     settings.use_fesd = 0;
    [results] = local_minima_experiment(scenario,settings,model);
    % Scenario-3  : use FESD with homotopy
    scenario.scenario_name = 'local_min_GL4_fesd_homotopy';
    settings.use_fesd = 1;
    settings.sigma_0 = 1e0;
    settings.sigma_N = 1e-15;
    settings.N_homotopy = 15;
    [results] = local_minima_experiment(scenario,settings,model);
    % Scenario-4  : use Standard with homotopy
    scenario.scenario_name = 'local_min_GL4_std_homotopy';
    settings.use_fesd = 0;
    [results] = local_minima_experiment(scenario,settings,model);
end
%% Objective function experimets
if run_objective_function_experiment
    % Scenario 1 : FESD
    scenario.scenario_name = 'objective_GL4_fesd_fixed';
    settings.use_fesd = 1;
    settings.mpcc_mode = 1;
    settings.sigma_0 = 1e-15;
    settings.sigma_N = 1e-15;
    settings.N_homotopy = 3;
    [results] = objective_function_experiment(scenario,settings,model);
    % Scenario 2 : Standard
    scenario.scenario_name = 'objective_GL4_std_fixed';
    settings.use_fesd = 0;
    settings.mpcc_mode = 3;
    settings.sigma_0 = 1e-15;
    settings.sigma_N = 1e-15;
    settings.N_homotopy = 1;
    [results] = objective_function_experiment(scenario,settings,model);
end