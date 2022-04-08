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
%% NOS-NOC settings
[settings] = default_settings_fesd();  %% Optionally call this function to have an overview of all options.
settings.time_freezing = 1; 
settings.n_s = 3; 
settings.print_level = 3;
%% model equations
model.T = 4; model.N_stages = 15; model.N_finite_elements = 3;
model.x0 = [0;0.5;0;0;0];
q = MX.sym('q',2); v = MX.sym('v',2); t = MX.sym('t');

model.x = [q;v;t];
model.c = q(2); model.S = [1; -1];
u = MX.sym('u',2); model.u = u;

f_1 = [v;u(1);u(2)-9.81;1]; f_2 = [0;v(2);0;-10*(q(2))-0.211989*v(2);0];
model.F = [f_1 f_2];
%% Objective and constraints
model.f_q = u'*u; model.f_q_T = 10*v'*v;
model.g_ineq = u'*u-7^2;
model.g_terminal = q-[4;0.25];
%% Solve and plot
[solver,solver_initalization, model,settings] = create_nlp_fesd(model,settings);
[results,stats] = homotopy_solver(solver,model,settings,solver_initalization);
plot_result_ball(model,settings,results,stats)
fprintf('Objective values is: %2.4f \n',full(results.f));
