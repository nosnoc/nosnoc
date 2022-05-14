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
clear all
clc
close all
import casadi.*

%% Settings
[settings] = default_settings_nosnoc();
settings.n_s = 2;                       
settings.mpcc_mode = 6; 
settings.s_elastic_max = 1e0;             
settings.opts_ipopt.ipopt.max_iter = 1e3;
% settings.N_homotopy = 11;
% Step Equlibration
settings.step_equilibration_penalty = 1;
settings.step_equilibration_mode = 3; %1 - penality relaxation, 3 - step equilibration as hard constraints
%% Time settings
settings.time_freezing = 1;
settings.time_optimal_problem = 1;
% Time freezing scaling / Speed of Time
settings.s_sot_max = 10;
settings.use_speed_of_time_variables =  1; % introduce s_tof for e
settings.local_speed_of_time_variable = 1;
% Grid settings
settings.equidistant_control_grid = 1;
settings.stagewise_clock_constraint = 1;
%% Model Settings
model.fuel_cost_on = 0;
% Discretization parameters
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