function [output] = solve_car_turbo_with_nosnoc(model,settings)
%
%    This file is part of NOSNOC.
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
%% settings
import casadi.*
%%
N_trails = settings.N_trails;
q_goal = 200;
v_goal = 0;
v_max = 25; % maximal velocity
u_max = 5;
v_trash_hold = 10;
%% Model - define all problem functions and
% Discretization parameters
% Symbolic variables and bounds
q = SX.sym('q'); % position
v = SX.sym('v'); % velocity
model.x = [q;v]; % add all important data to the struct model,
model.x0 = [0;0]; % inital value

model.lbx = [-inf;-v_max];
model.ubx = [inf;v_max];

% control
u = SX.sym('u');
model.u = u;
% bounds
model.lbu = -u_max;
model.ubu = u_max;
% Dyanmics and the regions
f_1 = [v;u]; % mode 1 - nominal
f_2 = [v;3*u]; % mode 2 - turbo
model.F = [f_1 f_2];
% Constraint for the regions
model.c = v-v_trash_hold;
model.S = [-1;1];

% Objective
if ~settings.time_optimal_problem
    % note that if q_goal is increased, model.T should be increased as
    % otherwise the problem might become infeasible
    % for time optimal problems we do not have to provide an objective
    model.f_q = u^2;
end
% terminal constraint
model.g_terminal = [q-q_goal;v-v_goal];
%% Solve OCP
% This functions formulates and discretized the OCP. We obtain an matheatmical programm with complementarity constraint which is solved  in a homotopy procedure.
cpu_time_all = [];
error_all = [];
for ii = 1:N_trails
[results,stats,model,settings] = nosnoc_solver(model,settings);
cpu_time_all = [cpu_time_all, stats.cpu_time_total];
results.T_opt = results.t_grid(end);
[tout,yout,error]= car_turbo_sim(results.u_opt,results.T_opt,model.N_stages,1);
error_all = [error_all ;error];
end
results.T_opt = results.t_grid(end);
[tout,yout,error]= car_turbo_sim(results.u_opt,results.T_opt,model.N_stages,1);


output.T_opt = results.T_opt;
output.error = error;
output.error_all= error_all;
output.tout = tout;
output.yout = yout;
output.results = results;
output.cpu_time_all = cpu_time_all;
output.cpu_time =  mean(cpu_time_all);
output.N_stages =  model.N_stages;
output.N_finite_elements =  model.N_finite_elements ;



end

