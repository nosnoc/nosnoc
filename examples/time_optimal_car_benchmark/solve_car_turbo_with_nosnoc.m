function [output] = solve_car_turbo_with_nosnoc(model,settings)
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

%% settings
import casadi.*
%%
N_trails = settings.N_trails;
q_goal = 200;
v_goal = 0;
v_max = 25; % maximal velocity
u_max = 5;
v_trash_hold = 10;
%% Model - define all problem functions and
% Discretization parameters
% Symbolic variables and bounds
q = SX.sym('q'); % position
v = SX.sym('v'); % velocity
model = NosnocModel();
model.x = [q;v]; % add all important data to the struct model,
model.x0 = [0;0]; % inital value

model.lbx = [-inf;-v_max];
model.ubx = [inf;v_max];

% control
u = SX.sym('u');
model.u = u;
% bounds
model.lbu = -u_max;
model.ubu = u_max;
% Dyanmics and the regions
f_1 = [v;u]; % mode 1 - nominal
f_2 = [v;3*u]; % mode 2 - turbo
model.F = [f_1 f_2];
% Constraint for the regions
model.c = v-v_trash_hold;
model.S = [-1;1];

% Objective
if ~settings.time_optimal_problem
    % note that if q_goal is increased, model.T should be increased as
    % otherwise the problem might become infeasible
    % for time optimal problems we do not have to provide an objective
    model.f_q = u^2;
end
% terminal constraint
model.g_terminal = [q-q_goal;v-v_goal];
%% Solve OCP
% This functions formulates and discretized the OCP. We obtain an matheatmical programm with complementarity constraint which is solved  in a homotopy procedure.
cpu_time_all = [];
error_all = [];
for ii = 1:N_trails
solver = NosnocSolver(model, settings);
[results,stats] = solver.solve();

cpu_time_all = [cpu_time_all, stats.cpu_time_total];
results.T_opt = results.t_grid(end);
[tout,yout,error]= car_turbo_sim(results.u_opt,results.T_opt,settings.N_stages,1);
error_all = [error_all ;error];
end
results.T_opt = results.t_grid(end);
[tout,yout,error]= car_turbo_sim(results.u_opt,results.T_opt,settings.N_stages,1);


output.T_opt = results.T_opt;
output.error = error;
output.error_all= error_all;
output.tout = tout;
output.yout = yout;
output.results = results;
output.cpu_time_all = cpu_time_all;
output.cpu_time =  mean(cpu_time_all);
output.N_stages =  settings.N_stages;
output.N_finite_elements =  settings.N_finite_elements ;



end

