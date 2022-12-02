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

%% Settings
[settings] = default_settings_nosnoc();
settings.n_s = 3;                       
settings.opts_ipopt.ipopt.max_iter = 1e3;
%% Time settings
settings.time_freezing = 1;
settings.time_freezing_hysteresis  = 1;
settings.time_optimal_problem = 1;
% Time-freezing scaling / speed of time
settings.s_sot_max = 10;
settings.s_sot_min = 1;
settings.rho_sot = 1e-1;
settings.use_speed_of_time_variables = 1; 
settings.local_speed_of_time_variable = 1;
% solver settings
settings.opts_ipopt.ipopt.tol = 1e-8;
settings.comp_tol = 1e-8;
settings.mpcc_mode = 3; 
settings.cross_comp_mode = 8;
settings.homotopy_update_rule = 'superlinear';

%% Model Settings
model.fuel_cost_on = 0;
model.N_finite_elements = 3;
model.N_stages = 10;
model.T = 5;
%% solve OCP
model = car_hystheresis_model_voronoi(model);
[results,stats,model,settings] = nosnoc_solver(model,settings);
%% Read and plot Result
plot_results_car_hysteresis(results,settings,model,stats)   
%%
T_opt = results.T_opt;
N_stages = model.N_stages;
u_opt = results.u_opt;
[tout,yout,error] = car_hysteresis_sim(u_opt,T_opt,N_stages);