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
%
clear all
clc 
close all
import casadi.*

%% Info
% This is an example from 
% Stewart, D.E., 1996. A numerical method for friction problems with multiple contacts. The ANZIAM Journal, 37(3), pp.288-308.
% It considers 3 independent switching functions and it demonstrates the
% generalization of the FESD scheme presented in the NOSNOC software parep
%% settings
% collocation settings
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();


problem_options.n_s = 2;                            
problem_options.rk_scheme = RKSchemes.RADAU_IIA;     

solver_options.print_level = 2;
solver_options.s_elastic_max = 1e1;                    
solver_options.sigma_0 = 1;
solver_options.complementarity_tol = 1e-9;
solver_options.N_homotopy = 10;
% solver_options.homotopy_update_rule = 'superlinear';
% solver_options.homotopy_update_slope = 0.2;
% solver_options.homotopy_update_exponent = 2.5;
%% Generate Model
model = blocks_with_friction();
%% Simulation settings
N_finite_elements = 3;
T_sim = 12;
N_sim = 85;

problem_options.dcs_mode = 'Stewart';
% problem_options.dcs_mode = 'Heaviside';

problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_finite_elements;
problem_options.N_sim = N_sim;
solver_options.use_previous_solution_as_initial_guess = 1;
%% Call FESD Integrator
integrator = nosnoc.Integrator(model, problem_options, solver_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();
%% Get variables into main workspace
if problem_options.dcs_mode == 'Stewart'
    theta_res = integrator.get("theta");
    lambda_res = integrator.get("lambda");
    mu_res = integrator.get("mu");
else
    alpha_res = integrator.get("alpha");
    lambda_n_res = integrator.get("lambda_n");
    lambda_p_res = integrator.get("lambda_p");
end
h_res = integrator.get("h");
plot_results_friction_blocks

