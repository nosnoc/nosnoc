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
clear all
clc
close all
import casadi.*
%%
switching_case = 'sliding_mode';
%  Options: 'crossing' 'sliding_mode', 'spontaneous_switch' , 'leave_sliding_mode',
%% NOSNOC settings
[settings] = default_settings_nosnoc();  %% Optionally call this function to have an overview of all options.
settings.n_s = 1;
settings.use_fesd = 1;
settings.mpcc_mode = 3;
settings.kappa = 0.1;
settings.irk_scheme = 'Gauss-Legendre';
settings.irk_scheme = 'Radau-IIA';
settings.irk_representation= 'differential';
settings.irk_representation= 'integral';
settings.pss_mode = 'Step';
settings.cross_comp_mode = 3;
settings.lift_irk_differential = 1;
settings.print_level = 5;
settings.step_equilibration = 0;
settings.heuristic_step_equilibration = 0;
% settings.gamma_h = 0.9;

% discretization parameters
N_sim = 1;
T_sim = 0.75;

model.N_sim = N_sim;
model.N_finite_elements = 2;
model.T_sim = T_sim;

model.x0 = -0.50;
x = SX.sym('x',1);
model.x = x;
model.c = x;
model.S = [-1; 1];
f_1 = [1]; f_2 = [-1];
model.F = [f_1 f_2];
[results,stats,model] = integrator_fesd(model,settings);
%
figure
plot(results.t_grid,results.x_res)
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$x(t)$','Interpreter','latex')
grid on



