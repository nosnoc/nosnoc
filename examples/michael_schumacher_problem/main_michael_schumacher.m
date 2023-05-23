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

% This file is part of NOSNOC.%
clear all; clc; close all
import casadi.*
%%
% An optimal control problem from: 
% Optimal control of systems with discontinuous differential equations 
% April 2012, Numerische Mathematik 114(4):653-695
% DOI: 10.1007/s00211-009-0262-2
% David E. Stewart Mihai Anitescu
%% A very difficult problem is with:
% path_constraint = 'track';
% settings.time_optimal_problem = 1;

%% NOSNOC settings
path_constraint = 'linear';
track_width = 0.5;

settings = NosnocOptions();  %% Optionally call this function to have an overview of all options.
settings.time_optimal_problem = 0;
settings.n_s = 2; 
settings.sigma_0 = 1e0;
settings.use_fesd = 1;
settings.cross_comp_mode = 3;
settings.T_final_max = 5*pi;
settings.T_final_min = 2;
settings.opts_casadi_nlp.ipopt.linear_solver = 'ma57';
settings.dcs_mode = 'Step';

% settings.relax_terminal_constraint = 2;
% settings.rho_terminal = 1e3;

%% model equations
model = NosnocModel();
model.x0 = zeros(5,1);
model.T = 5; 
model.N_stages = 25; 
model.N_finite_elements = 2;

%% Variable defintion
% Declare model variables
qx = SX.sym('qx'); qy = SX.sym('qy');
vx = SX.sym('vx'); vy = SX.sym('vy');
alpha = SX.sym('alpha');
tangent = [cos(alpha);sin(alpha)];
normal = [-sin(alpha);cos(alpha)];
q = [qx;qy];
v = [vx;vy];
x = [q;v;alpha];
model.x = x;
model.lbx = [0;0;-inf;-inf;-inf];

%%  control
a = SX.sym('a');
s = SX.sym('s');
u = [a;s];
u_max  = 2;
model.u = u;
model.lbu = -u_max*ones(2,1);
model.ubu = u_max*ones(2,1);
model.u0 = [u_max;0];
%% modes of the PSS
Friction_max = 4;
f_1 = [v;a*tangent+Friction_max*normal;s*(tangent'*v)];
f_2 = [v;a*tangent-Friction_max*normal;s*(tangent'*v)];
% Switching functions and modes
model.c = normal'*v;
model.S = [1;-1];
model.F = [f_1 f_2];
%%  general nonlinear constinrst
% model.g_path = qy-sin(qx);
switch path_constraint 
    case 'none'
%         model.g_path = [];
        q_target = [3*pi;0];

    case 'linear'
        model.g_path = qy-(qx);
        q_target = [3*pi;3*pi];

    case 'nonlinear'
        omega = 3/8;
        omega = 1; 
%         omega = 0.3;
        model.g_path = qy-sin(omega*qx);
        q_target = [3*pi;sin(omega*3*pi)];

    case 'track'
        arg1 = qx-pi;
        arg2 = qx-2*pi;
        sig = 1e-1;
        step1 = 0.5*(1+tanh(arg1/sig));
        step2 = 0.5*(1+tanh(arg2/sig));
        model.g_path = qy-(sin(qx)*(1-step1)+(pi-qx)*step1*(1-step2)+(-pi-sin(qx))*step2);
        q_target = [3*pi;-pi];
end
if ~isequal(path_constraint,'none')
    model.g_path_lb = [-track_width];
    model.g_path_ub = [+track_width];
end
%% objective
if ~settings.time_optimal_problem 
    f_q = u'*u;
    model.T = 5;
else
    model.f_q = 0.01*u'*u;
end

% Terminal Constraint
model.g_terminal = [q-q_target];

%% Solve
solver = NosnocSolver(model, settings);
[results,stats] = solver.solve();
fprintf('Objective values is: %2.4f \n',full(results.f_opt));
%% plot
plot_results_ms
