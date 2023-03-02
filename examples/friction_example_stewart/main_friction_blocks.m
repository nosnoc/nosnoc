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
[settings] = default_settings_nosnoc();  %% Optionally call this function to have an overview of all options.
settings.n_s = 2;                            
settings.irk_scheme = 'Radau-IIA';     
settings.irk_representation = 'differential';
settings.lift_irk_differential = 1;
settings.print_level = 2;
settings.use_fesd = 1;
settings.mpcc_mode = 'Scholtes_ineq';
settings.sigma_0 = 1e4;
settings.s_elastic_max = 1e1;                    
settings.cross_comp_mode = 3;
settings.comp_tol = 1e-9;
settings.N_homotopy = 10;
settings.homotopy_update_rule = 'superlinear';
settings.homotopy_update_slope = 0.2;
settings.homotopy_update_exponent = 2.5;
%% Generate Model
model = blocks_with_friction();
%% Simulation setings
N_finite_elements = 3;
T_sim = 12;
% T_sim = 4;
N_sim = 85;

settings.dcs_mode = 'Stewart';
% settings.dcs_mode = 'Step';

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

