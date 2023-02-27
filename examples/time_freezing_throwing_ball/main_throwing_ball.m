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
%
clear all
clc
close all
import casadi.*
%% model parameters
e = 0.9; u_max = 9; beta = 0.0; 
%% NOSNOC settings
[settings] = default_settings_nosnoc();  %% Optionally call this function to have an overview of all options.
settings.time_freezing = 1; 
settings.n_s = 3; 
settings.mpcc_mode = 'Scholtes_ineq';
settings.homotopy_update_rule = 'superlinear';
settings.step_equilibration = 'direct_homotopy';
q_target = [4;0.5];
model.T = 4; 
model.N_stages = 20; 
model.N_finite_elements = 3;
% model equations
q = SX.sym('q',2);
v = SX.sym('v',2); 
u = SX.sym('u',2);
model.x0 = [0;0.5;0;0];
model.x = [q;v]; model.u = u; model.e = e ;
model.f_c = q(2);
model.f_v = [u-[0;9.81]-beta*v*sqrt(v(1)^2^2+v(2)^2+1e-3)]; 
% Objective and constraints
model.f_q = u'*u; model.f_q_T = 100*v'*v;
model.g_path = u'*u-u_max^2;
model.g_terminal = q-[q_target];
[results,stats,model,settings] = nosnoc_solver(model,settings);
%%
plot_result_ball(model,settings,results,stats)
fprintf('Objective values is: %2.4f \n',full(results.f_opt));
fprintf('Final time is: %2.4f \n',full(results.T_opt));
if isempty(results.T_opt)
    results.T_opt = results.t_grid(end);
end
[tout,yout,error] = bouncing_ball_sim(results.u_opt,results.T_opt,model.N_stages,model.x0(1:4),beta,0.9,q_target);
