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
function [model] = car_hystheresis_model(model_in)

import casadi.*
unfold_struct(model_in,'caller');
%% Discretization parameters
N_stages = 30;
N_finite_elements = 3;

active_control = 1; % for ocps

smooth_model = 0;
time_freezing_in_model = time_freezing;
%% goal
q_goal = 200;
v_goal = 0;
%% Bounds
v_max = 20;
umax = 5;

%% Time freezing and hysteresis settings
% normal ocp
if time_optimal_problem 
    N_stages = 35;
    N_finite_elements = 3;


    N_stages = 20;
    N_finite_elements = 3;
    % time optimal OcP
    fuel_cost_on = 0;
    control_cost_on = 0;
    v_trash_hold = 10; % switching criterion
    T = 3;   
else
    N_stages = 20;
    N_finite_elements = 3;
    fuel_cost_on = 0;
    control_cost_on = 1;
    v_trash_hold = 10; % switching criterion
    T = 15; % (here determined latter depeding on omega)
end


%% Model Parameters


% inital value
q0 = 0;
v0 = 0;
L0 =  0; % cost integral
a0 = 0;
t0 = 0;


% fual costs of turbo and nominal
% case 1
Pn = 1;
Pt = 1;
% case 2
Pn = 10;
Pt = 1;
%% Inital Value
x0 = [q0;v0;L0;a0;t0];

%% Time horizon
h = T/N_stages;
%% Model parameters for time freezing

%% Variable defintion
q = MX.sym('q');
v = MX.sym('v');
L = MX.sym('L');
a = MX.sym('a');
t = MX.sym('t');

x = [q;v;L;a;t];
n_x = length(x);


if active_control
    lbx = -[inf;0;inf;inf;inf];
    ubx = [inf;v_max;inf;inf;inf];
else
    lbx = -inf*ones(n_x,1);
    ubx = inf*ones(n_x,1);
end
% every constraint funcion corresponds to a simplex (note that the c_i might be vector valued)
c = v-v_trash_hold;

% sign matrix for the modes
S = [-1;...
       1];


%% control
u = MX.sym('u');
n_u = 1;
u0 = [0];
if active_control
    lbu = -umax*ones(n_u,1);
    ubu = umax*ones(n_u,1);
else
    lbu = umax*ones(n_u,1);
    ubu = umax*ones(n_u,1);
end

%% modes of the ODEs layers (for all  i = 1,...,n_simplex);


if smooth_model
    f_nominal_1 = [v;3*u;Pn;0;1];
    f_nominal_2 = [v;3*u;Pn;0;1];
else
    f_nominal_1 = [v;u;Pn;1;1];
    f_nominal_2 = [v;3*u;Pt;0;1*(1-time_freezing_in_model)];
    f_nominal_2 = [v;3*u;Pt;0;1];
end

f_11 = f_nominal_1;
f_12 = f_nominal_2;

% in matrix form
f_1 = [f_11 f_12];
F = f_1;
%% objective

f_q = fuel_cost_on*L+control_cost_on*u^2;
% Terminal Cost
f_q_T = 0;

%%  general nonlinear constinrst
general_nonlinear_constraint  = 0;
g_ineq = u^2;
g_ineq_lb = [-inf];
g_ineq_ub = [inf];
g_ineq_fun  = Function('g_ineq_fun',{x,u},{g_ineq});

%% terminal constraint
terminal_constraint = 1;

terminal_constraint =terminal_constraint*active_control;
g_terminal = [q-q_goal;v-v_goal];
g_terminal_lb = [0;0];
g_terminal_ub = [0;0];
g_terminal_fun = Function('g_terminal_fun',{x},{g_terminal});

%% Generic part
% (make of local workspace a struct and pass to output
names = who;
for ii = 1:length(names)
    eval([ 'model.' names{ii} '=' names{ii} ';'])
end
end

