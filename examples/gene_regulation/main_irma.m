% BSD 2-Clause License

% Copyright (c) 2023, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

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

%% Info
% Example synthetic benchmark from:
% Numerical simulation of piecewise-linear models of gene regulatory networks using complementarity systems
% V. Acary, H. De Jong, B. Brogliato
%%
clear; clc; close all
import casadi.*
%% Model parameters
switch_on = 1;
lifting = true;

%% Discretization
N_finite_elements = 4;
T_sim = 1000;
N_sim = 2000;


T_sim = 1000/8;
N_sim = 2000/8;

%% Settings
problem_options = nosnoc.Options();
integrator_options = nosnoc.integrator.Options();
solver_options = integrator_options.fesd_solver_opts;
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.print_level = 2;
problem_options.n_s = 2;
solver_options.homotopy_update_rule = 'superlinear';
solver_options.complementarity_tol = 1e-5;
solver_options.homotopy_steering_strategy = 'ELL_INF';
solver_options.decreasing_s_elastic_upper_bound = true;

% Generate model
model = irma_model(switch_on, lifting);

% Discretization settings
problem_options.N_finite_elements = N_finite_elements;
problem_options.T_sim = T_sim;
problem_options.N_sim = N_sim;

integrator = nosnoc.Integrator(model, problem_options, integrator_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();

results.x = x_res;
results.t_grid = t_grid;
results.alpha = integrator.get('alpha');
plot_irma(results);
