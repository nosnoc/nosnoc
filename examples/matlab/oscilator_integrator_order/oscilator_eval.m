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
function [U,V] = oscilator_eval(X,Y)
%OSCILATOR_EVAL Summary of this function goes here
%   Detailed explanation goes here
omega = 2*pi;
A1 = [1 omega;...
    -omega 1];
A2 = [1 -omega;...
    omega 1];
U = [];
V = [];
for ii=1:size(X,1)
    x = X(ii,:);
    y = Y(ii,:);
    u = 0.5*(1+sign(x.^2+y.^2-1)).*(A1(1,1).*x+A1(1,2).*y) + 0.5*(1-sign(x.^2+y.^2-1)).*(A2(1,1).*x+A2(1,2).*y);
    v = 0.5*(1+sign(x.^2+y.^2-1)).*(A1(2,1).*x+A1(2,2).*y) + 0.5*(1-sign(x.^2+y.^2-1)).*(A2(2,1).*x+A2(2,2).*y);
    U = [U;u];
    V = [V;v];
end


end

