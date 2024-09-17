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
clear; clc; close all;
import casadi.*
%% Experiments for derivatives of objective and local minima
run_objective_function_experiment = 1;
run_local_minima_experiment = 0;

N_samples = 50;
scenario.N_samples = N_samples;
scenario.x0_vec = linspace(-2,-0.8,N_samples);
scenario.save_results = 0;
%% load default nosnoc options
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();

% modify options
problem_options.n_s = 2;
problem_options.rk_scheme = RKSchemes.GAUSS_LEGENDRE;
problem_options.equidistant_control_grid = 0;
problem_options.step_equilibration = StepEquilibrationMode.direct;
problem_options.T = 2;
problem_options.N_stages = 1;
problem_options.N_finite_elements = 25;

%% nosnoc model
model = nosnoc.model.Pss();
model.x0 = -1;
% Variable defintion
x = SX.sym('x');
model.x = x;
model.c = x;
model.S = [1;-1];
% modes of the ODE
f_1 = 1;
f_2 = 3;
model.F = [f_1 f_2];
% objective
model.f_q = x^2;
model.f_q_T = (x-5/3)^2;
if run_local_minima_experiment
    % Scenario-1  : use FESD with fixed parameters
    problem_options.use_fesd = 1;
    solver_options.sigma_0 = 1e-15;
    solver_options.sigma_N = 1e-15;
    solver_options.N_homotopy = 1;
    scenario.scenario_name = ['local_min_fesd_fixed'];
    [results] = local_minima_experiment(scenario, problem_options, solver_options, model);

    % Scenario-2  : use Standard with fixed parameters
    problem_options.use_fesd = 0;
    scenario.scenario_name = ['local_min_std_fixed'];
    [results] = local_minima_experiment(scenario, problem_options, solver_options, model);

    % Scenario-3  : use FESD with homotopy
    scenario.scenario_name = ['local_min_fesd_homotopy'];
    problem_options.use_fesd = 1;
    solver_options.sigma_0 = 1e0;
    solver_options.sigma_N = 1e-15;
    solver_options.N_homotopy = 15;
    [results] = local_minima_experiment(scenario, problem_options, solver_options, model);

    % Scenario-4  : use Standard with homotopy
    scenario.scenario_name = ['local_min_std_homotopy'];
    problem_options.use_fesd = 0;
    [results] = local_minima_experiment(scenario, problem_options, solver_options, model);
end
%% Objective function experimets
if run_objective_function_experiment
    % Scenario 1 : FESD
    scenario.scenario_name = 'objective_fesd_fixed';
    problem_options.use_fesd = 1;
    settings.mpcc_mode = 'direct';
    settings.sigma_0 = 1e-15;
    settings.sigma_N = 1e-15;
    settings.N_homotopy = 3;
    results = objective_function_experiment(scenario, problem_options, solver_options, model);

    % Scenario 2 : Standard
    scenario.scenario_name = 'objective_std_fixed';
    problem_options.use_fesd = 0;
    settings.mpcc_mode = 'scholtes_ineq';
    settings.sigma_0 = 1e-15;
    settings.sigma_N = 1e-15;
    settings.N_homotopy = 1;
    results = objective_function_experiment(scenario, problem_options, solver_options, model);
end
