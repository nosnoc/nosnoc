% BSD 2-Clause License

% Copyright (c) 2022, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% This file is part of NOSNOC.
%
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
settings = default_settings_nosnoc();
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