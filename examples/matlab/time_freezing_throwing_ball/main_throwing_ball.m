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
%% model parameters
e = 0.9; u_max = 9; beta = 0.0; 
%% NOSNOC settings
[settings] = default_settings_nosnoc();  %% Optionally call this function to have an overview of all options.
settings.time_freezing = 1; 
settings.n_s = 3; 
q_target = [4;0.5];
model.T = 4; 
model.N_stages = 20; 
model.N_finite_elements = 3;
% model equations
% q = MX.sym('q',2); v = MX.sym('v',2); u = MX.sym('u',2);
q = SX.sym('q',2); v = SX.sym('v',2); u = SX.sym('u',2);
model.x0 = [0;0.5;0;0];
model.x = [q;v]; model.u = u; model.e = e ;
model.c = q(2);
model.f = [u-[0;9.81]-beta*v*sqrt(v(1)^2^2+v(2)^2+1e-3)]; 
% Objective and constraints
model.f_q = u'*u; model.f_q_T = 100*v'*v;
model.g_ineq = u'*u-u_max^2;
model.g_terminal = q-[q_target];
[results,stats,model,settings] = nosnoc_solver(model,settings);
%%
plot_result_ball(model,settings,results,stats)
fprintf('Objective values is: %2.4f \n',full(results.f_opt));
fprintf('Final time is: %2.4f \n',full(results.T_opt));
if isempty(results.T_opt)
    results.T_opt = results.t_grid(end);
end
[tout,yout,error] = bouncing_ball_sim(results.u_opt,results.T_opt,model.N_stages,model.x0(1:4),beta,0.9,q_target);
