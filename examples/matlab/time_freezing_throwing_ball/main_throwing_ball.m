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
settings.n_s = 4; 
settings.print_level = 3;
% settings.irk_scheme = 'Gauss-Legendre';
% settings.irk_representation
settings.mpcc_mode = 3;
settings.use_fesd  = 1;
% settings.N_homotopy = 10;
settings.time_optimal_problem = 1;
u_max = 9;
x_target = 4;
model.x_target = x_target;
%% model equations
model.T = 2; 
model.N_stages = 5; 
model.N_finite_elements = 3;
model.x0 = [0;0.5;0;0;0];
% model.x0 = [nan;nan;0;0;0]
q = MX.sym('q',2); v = MX.sym('v',2); t = MX.sym('t');
model.x = [q;v;t];
model.c = q(2); model.S = [1; -1];
u = MX.sym('u',2); model.u = u;

drag = sqrt(v(1)^2^2+v(2)^2+1e-3);
beta = 0.1;
f_1 = [v;u-[0;9.81]-beta*v*drag;1]; 
f_2 = [0;v(2);0;-10*(q(2))-0.211989*v(2);0];
model.F = [f_1 f_2];
%% Objective and constraints
model.f_q = u'*u*0; model.f_q_T = 10*v'*v;
if 1
    model.g_ineq = u'*u-u_max^2;
else
    model.lbu= -u_max*ones(2,1);
    model.ubu= u_max*ones(2,1);
end
model.g_terminal = q-[x_target;1];
%% Solve and plot
[results,stats,model,settings] = nosnoc_solver(model,settings);
% plot_result_ball(model,settings,results,stats)
%%
fprintf('Objective values is: %2.4f \n',full(results.f_opt));
fprintf('Final time is: %2.4f \n',full(results.T_opt));
if isempty(results.T_opt)
    results.T_opt = results.t_grid(end);
end
[tout,yout,error] = bouncing_ball_sim(results.u_opt,results.T_opt,model.N_stages,model.x0(1:4),beta);

%%

