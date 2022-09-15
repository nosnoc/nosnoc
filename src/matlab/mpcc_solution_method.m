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
%    NO-NOC is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with NOS-NOC; if not, write to the Free Software Foundation,
%    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%

function [nlp_data] = mpcc_solution_method(nlp_data,ind_data,mpcc_mode)
% Input
unfold_struct(nlp_data,'caller');
unfold_struct(ind_data,'caller');


% ind_theta = ind_z;
% ind_lambda = inz_zl
%%
n_mpcc = length(ind_g_mpcc);
switch mpcc_mode
    case 1
    case 2
    case 3

%         lbw(ind_theta) = -inf;
%         lbw(ind_lambda) = -inf;
% %         g = [g;z(ind_theta)-sigma_p*ones()];
%         g = [g;z(ind_lambda)-sigma_p];
%         lbg = [lbg; ]
%           g(ind_g_mpcc) = g(ind_g_mpcc)-ones(n_mpcc,1)*sigma_p;
%           lbg(ind_g_mpcc) = 0;
%           ubg(ind_g_mpcc) = 0;
    case 4
    case 5
end

%% Output
    nlp_data.J = J; nlp_data.g = g; nlp_data.lbg = lbg; nlp_data.ubg = ubg; 
    nlp_data.w = w; nlp_data.w0= w0; nlp_data.lbw = lbw; nlp_data.ubw = ubw;
    nlp_data.sigma_p = sigma_p;
end

