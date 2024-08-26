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

function [model] = two_gene_model(x0, lifting)
% Generate model for two gene regulatory network with given initial conditions
    import casadi.*
    model = nosnoc.model.Heaviside();
    % Initial Value
    model.x0 = x0;
    % Variables
    % Concentrations
    x = SX.sym('x', 2);
    % alphas for general inclusions
    alpha = SX.sym('alpha', 4);
    % Thresholds
    theta_1 = [4;8];
    theta_2 = [4;8];
    % Synthesis
    kappa = [40;40];
    % Degradation
    gamma = [4.5;1.5];
    % Switching function
    c = [x(1)-theta_1; x(2)-theta_2];
    % Switching multipliers
    s = [(1-alpha(2))*alpha(3);
         alpha(1)*(1-alpha(4))];
    if lifting
        beta = SX.sym('beta', 2);
        model.f_alg = beta - s; % lifted equations
        model.z = beta;

        % f_x
        f_x = -gamma.*x + kappa.*beta;
    else
        % f_x
        f_x = -gamma.*x + kappa.*s;
    end
    % set model parameters
    model.x = x;
    model.alpha = alpha;
    model.f_x = f_x;
    model.c = c;
end

