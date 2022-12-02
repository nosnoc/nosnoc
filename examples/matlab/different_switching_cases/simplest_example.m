% Copyright 2022 Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% This file is part of NOSNOC.

% The 2-Clause BSD License

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
clear all
clc
close all
import casadi.*
%%
switching_case = 'sliding_mode';
%  Options: 'crossing' 'sliding_mode', 'spontaneous_switch' , 'leave_sliding_mode',
%% NOSNOC settings
[settings] = default_settings_nosnoc();  %% Optionally call this function to have an overview of all options.
settings.n_s = 1;
settings.use_fesd = 1;
settings.mpcc_mode = 3;
settings.homotopy_update_slope = 0.1;
settings.irk_scheme = 'Gauss-Legendre';
settings.irk_scheme = 'Radau-IIA';
settings.irk_representation= 'differential';
settings.irk_representation= 'integral';
settings.pss_mode = 'Step';
settings.cross_comp_mode = 3;
settings.lift_irk_differential = 1;
settings.print_level = 5;
settings.step_equilibration = 0;
settings.heuristic_step_equilibration = 0;
% settings.gamma_h = 0.9;

% discretization parameters
N_sim = 1;
T_sim = 0.75;

model.N_sim = N_sim;
model.N_finite_elements = 2;
model.T_sim = T_sim;

model.x0 = -0.50;
x = SX.sym('x',1);
model.x = x;
model.c = x;
model.S = [-1; 1];
f_1 = [1]; f_2 = [-1];
model.F = [f_1 f_2];
[results,stats,model] = integrator_fesd(model,settings);
%
figure
plot(results.t_grid,results.x_res)
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$x(t)$','Interpreter','latex')
grid on



