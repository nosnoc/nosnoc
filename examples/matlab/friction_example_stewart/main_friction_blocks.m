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

%% Info
% This is an example from 
% Stewart, D.E., 1996. A numerical method for friction problems with multiple contacts. The ANZIAM Journal, 37(3), pp.288-308.
% It considers 3 independent switching functions and it demonstrates the
% generalization of the FESD scheme presented in the NOS NOC software parep
%% settings
% collocation settings
[settings] = default_settings_fesd();  %% Optionally call this function to have an overview of all options.
settings.kappa = 0.05;
settings.n_s = 2;                            % Degree of interpolating polynomial
settings.irk_scheme = 'radau';     % Collocation scheme: radau or legendre
settings.mpcc_mode = 3;                    % 1 - extact, 2 - smooth  ,3 -relax , 4 - penalty, 5 - l1 penalty
settings.s_elastic_max = 1e1;                    
settings.cross_comp_mode = 8;
%% Generate Model
model = blocks_with_friction();
%% Simulation setings
model.T_sim = 12;
model.N_stages = 2;
model.N_finite_elements = 1;
model.h = 0.05;
settings.gamma_h = 1;
model.T = model.N_stages*model.h;
settings.use_previous_solution_as_initial_guess = 1;
%% Call FESD Integrator
[results,stats,model] = integrator_fesd(model,settings);
%% Formulate NLP; 
%% Get variables into main workspace
unfold_struct(model,'base');
unfold_struct(settings,'base');
unfold_struct(results,'base');
plot_results_friction_blocks