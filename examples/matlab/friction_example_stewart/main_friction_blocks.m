%
%    This file is part of NOSNOC.
%
%    NOSNOC -- A software for NOnSmooth Numerical Optimal Control.
%    Copyright (C) 2022 Armin Nurkanovic, Moritz Diehl (ALU Freiburg).
%
%    NOSNOC is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 3 of the License, or (at your option) any later version.
%
%    NOSNOC is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with NOSNOC; if not, write to the Free Software Foundation,
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
% generalization of the FESD scheme presented in the NOSNOC software parep
%% settings
% collocation settings
[settings] = default_settings_nosnoc();  %% Optionally call this function to have an overview of all options.
settings.n_s = 2;                            
settings.irk_scheme = 'Radau-IIA';     
% settings.irk_scheme = 'Lobatto-IIIA';  
settings.irk_representation = 'differential';
settings.lift_irk_differential = 1;
settings.print_level = 2;
% settings.mpcc_mode = 5; % very fast    
settings.mpcc_mode = 3; % very robust e.g. with Gauss-Legendre for this example
settings.s_elastic_max = 1e1;                    
settings.cross_comp_mode = 3;
settings.comp_tol = 1e-6;
settings.N_homotopy = 7;
settings.homotopy_update_rule = 'superlinear';
settings.homotopy_update_slope = 0.2;
settings.homotopy_update_exponent = 2;
%% Generate Model
model = blocks_with_friction();
%% Simulation setings
N_finite_elements = 3;
% T_sim = 12;
T_sim = 4;
N_sim = 85;
% N_sim = 3;
% 
settings.pss_mode = 'Stewart';
% settings.pss_mode = 'Step';

model.T_sim = T_sim;
model.N_finite_elements = N_finite_elements;
model.N_sim = N_sim;

settings.use_previous_solution_as_initial_guess = 1;
%% Call FESD Integrator
[results,stats,model] = integrator_fesd(model,settings);
%% Get variables into main workspace
unfold_struct(model,'base');
unfold_struct(settings,'base');
unfold_struct(results,'base');
plot_results_friction_blocks

