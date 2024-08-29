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

%% About
% This script provids a tutoral how to formulate a time-optimal control
% problme of a simple car with two modes of operation.
% When the car exceeds a certain velocity treshold, its acceleration
% improves by the factor of 3 (turbo mode). For more information see the
% github page, the manual and the NOS-NOC related publications.
%% settings
clear all; clc; close all; 
import casadi.*
[settings] = NosnocOptions();  % Optionally call this function to have an overview of all options. Missing settings are anyway filled in latter with their respecitve values.
%% Choosing the Runge - Kutta Method and number of stages
settings.print_level = 3;
settings.n_s = 2;
settings.cross_comp_mode = 3;
settings.mpcc_mode = 'elastic_ineq';
settings.use_fesd = 1;
settings.use_speed_of_time_variables = 1;
settings.homotopy_update_rule = 'superlinear';

%% Time settings
settings.time_optimal_problem = 1;
%% Model - define all problem functions and
% Discretization parameters
settings.N_stages = 10; % number of control intervals
settings.N_finite_elements = 3; % number of finite element on every control interval (optionally a vector might be passed)
problem_options.T = 1;    % Time horizon

% Symbolic variables and bounds
q = SX.sym('q'); % position
v = SX.sym('v'); % velocity
model.x = [q;v]; % add all important data to the struct model,
model.x0 = [0;0]; % inital value

v_max = 20; % maximal velocity
model.lbx = [-inf;-v_max];
model.ubx = [inf;v_max];

% control
u = SX.sym('u');
model.u = u;
% bounds
u_max = 5;
model.lbu = -u_max;
model.ubu = u_max;

% Dyanmics and the regions
f_1 = [v;u]; % mode 1 - nominal
f_2 = [v;3*u]; % mode 2 - turbo
model.F = [f_1 f_2];
% Constraint for the regions
v_trash_hold = 10;
model.c = v-v_trash_hold;
model.S = [-1;1];

% Objective
if ~settings.time_optimal_problem
    % note that if q_goal is increased, problem_options.T should be increased as
    % otherwise the problem might become infeasible
    % for time optimal problems we do not have to provide an objective
    model.f_q = u^2;
end

%%  general nonlinear constraints are possible as well
% model.g_path = u^2;
% model.g_path_lb = [-inf];
% model.g_path_ub = [u_max^2];

%  model.g_path = v;
% model.g_path_lb = [-inf];
% model.g_path_ub = [20];
% settings.g_path_at_fe = 1;
% settings.g_path_at_stg = 1;
%% terminal constraint
q_goal = 200;
v_goal = 0;
% Add terminal constraint, if no upper and lower bound are provided, they are set to zero
model.g_terminal = [q-q_goal;v-v_goal];

%% Solve OCP
% This functions formulates and discretized the OCP. We obtain an matheatmical programm with complementarity constraint which is solved  in a homotopy procedure.
solver = NosnocSolver(model, settings);
[results,stats] = solver.solve();
plot_results_nosnoc_tutorial
if results.T_opt<1
   results.T_opt = results.t_grid(end) ;
end
 [tout,yout,error]= car_turbo_sim(results.u_opt,results.T_opt,settings.N_stages,1);
