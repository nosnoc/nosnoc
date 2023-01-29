% BSD 2-Clause License

% Copyright (c) 2022, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% This file is part of NOSNOC.

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

