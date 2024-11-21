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

clear all
clc
close all
import casadi.*


%% problem_options
% collocation problem_options
problem_options = nosnoc.Options();
problem_options.n_s = 2;                            % Degree of interpolating polynomial
problem_options.print_level = 2;
problem_options.T_sim = 4;
problem_options.N_sim = 40;
problem_options.N_finite_elements = 2;
%% Generate Model
model = temp_control_model_voronoi();
%% - solver_options settings
integrator_options = nosnoc.integrator.Options();
solver_options = integrator_options.fesd_solver_opts;
solver_options.homotopy_steering_strategy = 'ELL_INF';
solver_options.decreasing_s_elastic_upper_bound = true;
solver_options.N_homotopy = 6;
solver_options.print_level = 3;
solver_options.opts_casadi_nlp.ipopt.max_iter = 1e3;
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
solver_options.use_previous_solution_as_initial_guess = 1;
%% Call FESD Integrator
integrator = nosnoc.Integrator(model, problem_options, integrator_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();
%% Read and plot result
plot_results_for_paper
