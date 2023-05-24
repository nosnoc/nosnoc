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

function [model] = irma_model(switch_on, lifting)
% Generate model for two gene regulatory network with given initial conditions
    import casadi.*
    model = NosnocModel(); 
    % Initial Value
    model.x0 = [0.011;0.09;0.04;0.05;0.015];
    % Variables
    % Concentrations
    x = SX.sym('x', 5);
    % alphas for general inclusions
    alpha = SX.sym('alpha', 7);
    % Thresholds
    theta = {[0.01],[0.01;0.06;0.08],[0.035],[0.04],[0.01]};
    % Synthesis
    kappa = [1.1e-4, 9e-4;
             3e-4, 0.15;
             6e-4, 0.018;
             5e-4, 0.03;
             7.5e-4, 0.015];
    % Degradation
    gamma = [0.05;0.04;0.05;0.02;0.6];
    % Switching function
    c = [x(1)-theta{1};
         x(2)-theta{2};
         x(3)-theta{3};
         x(4)-theta{4};
         x(5)-theta{5}];
    if lifting
        if switch_on
            beta = SX.sym('beta', 1);
            model.g_z = beta - alpha(2)*(1-alpha(5));% lifted equations
            s = [1, alpha(6);
             1, alpha(1);
             1, alpha(3);
             alpha(2), beta(1);
             1, alpha(4)];
        else
            beta = SX.sym('beta', 2);
            model.g_z = beta - [alpha(1)*(1-alpha(7)); alpha(2)*(1-alpha(5))];% lifted equations
            s = [1, alpha(6);
             1, beta(1);
             1, alpha(3);
             alpha(2), beta(2);
             1, alpha(4)];
        end
        model.z = beta;
        f_x = -gamma.*x + sum(kappa.*s, 2);
    else
        % f_x
        s = [1, alpha(6);
             1, alpha(1)*(1-(1-switch_on)*(alpha(7)));
             1, alpha(3);
             alpha(2),alpha(2)*(1-alpha(5));
             1, alpha(4)];
        f_x = -gamma.*x + sum(kappa.*s, 2);
    end
    % set model parameters
    model.x = x;
    model.alpha = {alpha};
    model.f_x = f_x;
    model.c = c;
end

