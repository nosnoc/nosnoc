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
%% About
% This script provids a tutoral how to formulate a time-optimal control
% problme of a simple car with two modes of operation.
% When the car exceeds a certain velocity treshold, its acceleration
% improves by the factor of 3 (turbo mode). For more information see the
% github page, the manual and the NOS-NOC related publications.
%% settings
clear all; clc; close all; 
import casadi.*
[settings] = default_settings_fesd();  % Optionally call this function to have an overview of all options. Missing settings are anyway filled in latter with their respecitve values.
%% Choosing the Runge - Kutta Method and number of stages
settings.print_level = 3;
settings.irk_scheme = 'Gauss-Legendre';
settings.irk_scheme = 'Radau-IIA';
% settings.irk_scheme = 'Explicit-RK';
% settings.irk_scheme = 'Lobatto-IIIA';
settings.n_s = 2;
settings.mpcc_mode = 5;
settings.cross_comp_mode = 3;
% settings.s_elastic_0 = 1e1;
% settings.s_elastic_max = 1e4;
settings.pss_mode = 'Stewart';
% settings.pss_mode = 'Step';
settings.use_fesd = 1;

%% Time settings
% Here we can indicate tha the Optimal Control Problem (OCP) is a time optimal control problem so the
% solver can intorduce the needed time transfomrations and create the objective function.
settings.time_optimal_problem = 1;
%% Model - define all problem functions and
% Discretization parameters
model.N_stages = 10; % number of control intervals
model.N_finite_elements = 3; % number of finite element on every control intevral (optionally a vector might be passed)
model.T = 15;    % Time horizon

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
    % note that if q_goal is increased, model.T should be increased as
    % well, otherwise the problem might become infeasible
    %  for time optimal problems we do not have to provide an objective
    model.f_q = u^2;
end

%%  general nonlinear constraints are possible as well
% model.g_ineq = u^2;
% model.g_ineq_lb = [-inf];
% model.g_ineq_ub = [u_max^2];

%  model.g_ineq = v;
% model.g_ineq_lb = [-inf];
% model.g_ineq_ub = [20];
% settings.g_ineq_at_fe = 1;
% settings.g_ineq_at_stg = 1;
%% terminal constraint
q_goal = 200;
v_goal = 0;
% Add terminal constraint, if no upper and lower bound are provided, they are set to zero
model.g_terminal = [q-q_goal;v-v_goal];

%% Solve OCP
% This functions formulates and discretized the OCP. We obtain an matheatmical programm with complementarity constraint which is solved  in a homotopy procedure.
[results,stats,model,settings] = nosnoc_solver(model,settings);
plot_results_nosnoc_tutorial
