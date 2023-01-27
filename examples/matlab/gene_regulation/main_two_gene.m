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

%% Info
% TODO Add paper reference
%%
clear all
clc
close all
import casadi.*

%% Discretization
N_finite_elements = 2;
T_sim = 1;
N_sim = 10;

%% Settings
settings = default_settings_nosnoc();
settings.use_fesd = 1;
settings.irk_scheme = 'Radau-IIA';
settings.print_level = 2;
settings.n_s = 4;
settings.pss_mode = 'Step'; % General inclusions only possible in step mode.
settings.comp_tol = 1e-5;
settings.homotopy_update_rule = 'superlinear';

%% Generate different trajectories
results = [];
for x1 = 3:3:12
    for x2 = 3:3:12
        x0 = [x1;x2];
        % Generate model
        model = two_gene_model(x0);
        % Time
        model.N_finite_elements = N_finite_elements;
        model.T_sim = T_sim;
        model.N_sim = N_sim;

        [result,stats,model] = integrator_fesd(model,settings);
        results = [results,result];
    end
end

plot_two_gene(results, false)
