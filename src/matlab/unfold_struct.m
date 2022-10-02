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
%    License along with NOS-NOC; if not, write to the Free Software Foundation,
%    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%
function [varargout] = unfold_struct(input,mode)
% mode = 'caller'; % copy fields locally in the function where unfold_struct was called
% mode = 'base'; % copy fields  into the main workspace

%% 
import casadi.*
names = fieldnames(input);
for iii=1:length(names)
    eval([names{iii} '=input.' names{iii} ';']);
%     if eval(['~exist( '''  names{i} ''', ''var'')'])
        %         Check if a variable exists in the workspace, within a function
        assignin(mode, names{iii}, eval(names{iii}));
%     end
end
end




% for i=1:length(names)
%     eval([names{i} '=input.' names{i} ]);
%     eval(['varargout{' num2str(i) '} =' names{i}])
% end
