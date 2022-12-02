% Copyright 2022 Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% This file is part of NOSNOC.

% The 2-Clause BSD License

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
function [model] = temp_control_model_voronoi()
import casadi.*
%% Discretization parameters
N_stages = 2;
N_finite_elements = 1;
T = 0.1; % (here determined latter depeding on omega)
h = T/N_stages;
%% Model Parameters
active_control = 0;
% inital value
t0 = 0;
w0 = 0;
y0 = 15;

lambda_cool_down = -0.2; % cool down time constant of lin dynamics
u_heat = 10; % heater power

% jump points in x in the hysteresis function
y1 = 18;
y2 = 20;


z1 = [1/4;-1/4];
z2 = [1/4;1/4];
z3 = [3/4;3/4];
z4 = [3/4;5/4];
% Z = [1/4 1/4 3/4 3/4;...
%      -1/4 1/4 3/4 5/4]
Z = [z1 z2 z3 z4];

%% Inital Value
x0 = [y0;w0;t0];

%% Define model dimensions, equations, constraint functions, regions an so on.
n_simplex = 1;% number of Cartesian products in the model ("independet switches"), we call this layer
m_1 = 4;
m_vec = [m_1];
%% Variable defintion
y = MX.sym('y');
w = MX.sym('w');
t = MX.sym('t');

x = [y;w;t];
n_x = length(x);
lbx = -inf*ones(n_x,1);
ubx = inf*ones(n_x,1);

% linear transformation for rescaling of the switching function.
psi = (y-y1)/(y2-y1);
z = [psi;w];
% discriminant functions via voronoi

g_11 = norm([psi;w]-z1)^2;
g_12 = norm([psi;w]-z2)^2;
g_13 = norm([psi;w]-z3)^2;
g_14 = norm([psi;w]-z4)^2;



g_ind = [g_11;g_12;g_13;g_14];
c = g_ind;

%% control
u = MX.sym('u');
n_u = 1;
u0 = [0];

if active_control
    umax = 1e-3;
    lbu = -umax*ones(n_u,1);
    ubu = umax*ones(n_u,1);
else
    lbu = 0*ones(n_u,1);
    ubu = 0*ones(n_u,1);
end

%% modes of the ODEs layers (for all  i = 1,...,n_simplex);

u_heat = 10;
f_A = [lambda_cool_down*y+u_heat;0;1];
f_B = [lambda_cool_down*y;0;1];

a_push = 5;
f_push_down = [0;-a_push*(psi-1)^2/(1+(psi-1)^2);0];
f_push_up = [0;a_push*(psi)^2/(1+(psi)^2);0];

f_11 = 2*f_A-f_push_down;
f_12 = f_push_down;
f_13 = f_push_up;
f_14 = 2*f_B-f_push_up;
f_1 = [f_11 f_12 f_13 f_14];
F = f_1;
%% objective
% f_q = active_control*(u^2)+y^2;
%% Generic part
% (make of local workspace a struct and pass to output
names = who;
for ii = 1:length(names)
    eval([ 'model.' names{ii} '=' names{ii} ';'])
end
end

