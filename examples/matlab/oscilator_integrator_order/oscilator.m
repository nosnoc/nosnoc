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
function [model] = oscilator(model_in)

import casadi.*
%% Time horizon
if ~isempty(model_in)
    unfold_struct(model_in,'caller');
else
    disp('empty struct input');
    smooth_model = 0;
end
%% Model parameters

omega = 2*pi;

% omega = pi;
A1 = [1 omega;...
    -omega 1];
A2 = [1 -omega;...
    omega 1];
if smooth_model
    A2 = A1;
end
%% Inital Value
x0 = [exp(-1);0];
% x0 = [2*exp(-1);0];
if ~exist('R_osc')
    R_osc = 1;
end

%% Variable defintion
x1 = MX.sym('x1');
x2 = MX.sym('x2');
x = [x1;x2];
% every constraint funcion corresponds to a simplex (note that the c_i might be vector valued)
c = x1^2+x2^2-R_osc^2;
% sign matrix for the modes
S = [1;-1];
c = [c];

f_11 = A1*x;
f_12 = A2*x;
% in matrix form
F = [f_11 f_12];
%% Generic part
% (make of local workspace a struct and pass to output
names = who;

for ii = 1:length(names)
    eval([ 'model.' names{ii} '=' names{ii} ';'])
end

end

