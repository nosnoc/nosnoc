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
function [model] = blocks_with_friction()
import casadi.*


%% Initial value
x0 = [-1;1;-1;-1;1;1;0];
u0 = 0; % guess for control variables
%% Numer of ODE layers
% n_simplex = 3;% number of Cartesian products in the model ("independet switches"), we call this layer
% % number of modes in every simplex
% m_1 = 2;
% m_2 = 2;
% m_3 = 2;
% m_vec = [m_1 m_2 m_3];

%% Variable defintion
% differential states
q1 = SX.sym('q1');
q2 = SX.sym('q2');
q3 = SX.sym('q3');
v1 = SX.sym('v1');
v2 = SX.sym('v2');
v3 = SX.sym('v3');
t = SX.sym('t');

q = [q1;q2;q3];
v = [v1;v2;v3];
x = [q;v;t];

%% Control
% u = SX.sym('u');
% n_u = 1;  % number of parameters,  we model it as control variables and merge them with simple equality constraints
% 
% % Guess and Bounds
% u0 = 0;
% lbu  = -20*0;
% ubu  = 20*0;

%% Switching Functions
% every constraint function corresponds to a simplex (note that the c_i might be vector valued)
c1 = v1;
c2 = v2;
c3 = v3;
% sign matrix for the modes
S1 = [1;-1];
S2 = [1;-1];
S3 = [1;-1];
% discrimnant functions
S = {S1,S2,S3};
c = {c1,c2,c3};


%% Modes of the ODEs layers (for all  i = 1,...,n_simplex);
% part independet of the nonsmoothness
F_external = 0; % external force, e.g., control
F_input = 10; % variable force exicting
f_base = [v1;...
    v2;...
    v3;...
    (-q1)+(q2-q1)-v1;...
    (q1-q2)+(q3-q2)-v2;...
    (q2-q3)-v3+F_input*cos(pi*t);...
    ]/6;

f_base = [v1;...
    v2;...
    v3;...
    (-q1)+(q2-q1)-v1;...
    (q1-q2)+(q3-q2)-v2;...
    (q2-q3)-v3+F_external+F_input*(1*0+1*cos(pi*t));...
    1];
%
% for c1, h1,
f_11 = f_base+[0;0;0;-0.3;0;0;0];
f_12 = f_base+[0;0;0;+0.3;0;0;0];
% for c2, h2
f_21 = [0;0;0;0;-0.3;0;0];
f_22 = [0;0;0;0;0.3;0;0];
% for c3, h3
f_31 = [0;0;0;0;0;-0.3;0];
f_32 = [0;0;0;0;0;0.3;0];
% unfold_struct(model,'base');
% in matrix form
F1 = [f_11 f_12];
F2 = [f_21 f_22];
F3 = [f_31 f_32];

F = {F1 F2 F3};

%% Objective
% f_q = 0*u^2 + 1*v'*v;
% f_q_T = 10*v'*v;
%% Generic part
% (make of local workspace a struct and pass to output
names = who;
for ii = 1:length(names)
    eval([ 'model.' names{ii} '=' names{ii} ';'])
end

end