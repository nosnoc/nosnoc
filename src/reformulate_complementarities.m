% BSD 2-Clause License

% Copyright (c) 2023, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

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
function [g_out, lbg, ubg, cost] = reformulate_complementarities(g, mpcc_mode, sigma_p, s_elastic)
% For An explanation of the reformulations see the MpccMode class. 
% TODO mpcc_mode should be a structure with C-fun, relaxation type, etc.

    g_out = [];
    lbg = [];
    ubg = [];
    cost = 0;
    n_g = size(g,1);
    if mpcc_mode == MpccMode.Scholtes_ineq
        g_out = g - sigma_p;
        ubg = zeros(n_g,1);
        lbg = -inf * ones(n_g,1);
    elseif mpcc_mode == MpccMode.Scholtes_eq
        g_out = g - sigma_p;
        ubg = zeros(n_g,1);
        lbg = zeros(n_g,1);
    elseif mpcc_mode == MpccMode.ell_1_penalty
        cost = sum(g);
    elseif mpcc_mode == MpccMode.elastic_ineq
        g_out = g - s_elastic*ones(n_g,1);
        ubg = zeros(n_g,1);
        lbg = -inf * ones(n_g,1);
    elseif mpcc_mode == MpccMode.elastic_eq
        g_out = g - s_elastic*ones(n_g,1);
        ubg = zeros(n_g,1);
        lbg = zeros(n_g,1);
    elseif mpcc_mode == MpccMode.elastic_two_sided
        g_out = [g-s_elastic*ones(n_g,1);g+s_elastic*ones(n_g,1)];
        ubg = [zeros(n_g,1); inf*ones(n_g,1)];
        lbg = [-inf*ones(n_g,1);  zeros(n_g,1)];
    elseif mpcc_mode == MpccMode.elastic_ell_1_ineq
        g_out = g - s_elastic;
        ubg = zeros(n_g,1);
        lbg = -inf * ones(n_g,1);
    elseif mpcc_mode == MpccMode.elastic_ell_1_eq
        g_out = g - s_elastic;
        ubg = zeros(n_g,1);
        lbg = zeros(n_g,1);
    elseif mpcc_mode == MpccMode.elastic_ell_1_two_sided
        g_out = [g-s_elastic;g+s_elastic];
        ubg = [zeros(n_g,1); inf*ones(n_g,1)];
        lbg = [-inf*ones(n_g,1);  zeros(n_g,1)];
    end
end
