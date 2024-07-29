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

clear all
clc
close all
import casadi.*

%%
fuel_cost_on = 0;
%% Settings
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();

%% Options
problem_options.time_freezing = 1;
problem_options.time_freezing_hysteresis = 1;
%problem_options.time_optimal_problem = 1;
% Time-freezing scaling / speed of time
problem_options.s_sot_max = 100;
problem_options.s_sot_min = 0.1;
problem_options.rho_sot = 0;
problem_options.use_speed_of_time_variables = 1; 
problem_options.local_speed_of_time_variable = 1;
problem_options.stagewise_clock_constraint = 1;
problem_options.relax_terminal_constraint = ConstraintRelaxationMode.ELL_2;
problem_options.rho_terminal = 1e5;
% problem_options.relax_terminal_physical_time = ConstraintRelaxationMode.ELL_2;
% problem_options.rho_terminal_physical_time = 1e4;
% problem_options.relax_terminal_numerical_time = ConstraintRelaxationMode.ELL_1;
% problem_options.rho_terminal_numerical_time = 1e4;
problem_options.step_equilibration = StepEquilibrationMode.direct;
problem_options.rho_h = 1;
problem_options.n_s = 2;
problem_options.N_finite_elements = 3;
problem_options.N_stages = 10;
problem_options.T = 10;
problem_options.cross_comp_mode = 3;

% solver settings
solver_options.complementarity_tol = 1e-6;
solver_options.opts_casadi_nlp.ipopt.max_iter = 1e4;
solver_options.opts_casadi_nlp.ipopt.tol = 1e-7;
solver_options.opts_casadi_nlp.ipopt.acceptable_tol = 1e-5;
solver_options.opts_casadi_nlp.ipopt.acceptable_iter = 3;
solver_options.sigma_0 = 1;
solver_options.N_homotopy = 7;
solver_options.print_level = 3;
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';

%% Model
model = nosnoc.model.Pss();
%% Terminal constraint and bounds
q_goal = 150;
v_goal = 0;
v_max = 25;
u_max = 5;
%% Model Parameters
% Hystheresis parameters
v1 = 10;
v2 = 15;
% fual costs of turbo and nominal
Pn = 0;
Pt = 0;

%% Variable defintion
% states
q = SX.sym('q');
v = SX.sym('v');
L = SX.sym('L');
w = SX.sym('w');
t = SX.sym('t');
x = [q;v;L;w;t];
% controls
u= SX.sym('u');

model.x = x;
model.u = u;
%model.u0 = 5;

% Bounds on x and u
model.lbx = -[inf;v_max;inf;inf;inf];
model.ubx = [inf;v_max;inf;inf;inf];
model.lbu = -u_max;
model.ubu = u_max;
%% Inital Value
model.x0 = zeros(5,1);
% u0 = 10;
%% PSS via Voronoi Cuts
z1 = [1/4;-1/4];
z2 = [1/4;1/4];
z3 = [3/4;3/4];
z4 = [3/4;5/4];

psi = (v-v1)/(v2-v1);

g_1 = norm([psi;w]-z1)^2;
g_2 = norm([psi;w]-z2)^2;
g_3 = norm([psi;w]-z3)^2;
g_4 = norm([psi;w]-z4)^2;

model.g_ind = [g_1;g_2;g_3;g_4];

% modes of the ODEs layers
f_A = [v;u;Pn;0;1];
f_B = [v;3*u;Pt;0;1];

a_push = 1;
gamma_p = a_push*(psi^2/(1+psi^2));
gamma_n = a_push*((psi-1)^2/(1+(psi-1)^2));

f_push_down = [0;0;0;-gamma_n;0];
f_push_up = [0;0;0;gamma_p;0];

f_2 = f_push_down;
f_3 = f_push_up;
f_4 = 2*f_B-f_3;
f_1 = 2*f_A-f_2;
% in matrix form
model.F = [f_1 f_2 f_3 f_4];
%% objective and terminal constraint
model.f_q = fuel_cost_on*L;
% terminal constraint
model.g_terminal = [q-q_goal;v-v_goal];
%model.a_n = 1000;
% g_terminal_lb = zeros(2,1);
% g_terminal_ub = zeros(2,1);
% f_q_T = 1e3*[q-q_goal;v-v_goal]'*[q-q_goal;v-v_goal];
%% solve OCP
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();
%% Read and plot Result
x_res_full = ocp_solver.get_full("x");
theta_res_full = ocp_solver.get_full("theta");
lambda_res_full = ocp_solver.get_full("lambda");
h_res = ocp_solver.get("h");
u_res = ocp_solver.get("u")
t_res = x_res(end,:);
t_num = ocp_solver.get_time_grid();
t_num_full = ocp_solver.get_time_grid_full();


figure
subplot(4,1,1)
hold on
plot(t_num_full, theta_res_full(1,:))
plot(t_num_full, lambda_res_full(1,:))
for ii=1:length(t_num)
    xline(t_num(ii),'k:')
end
hold off
subplot(4,1,2)
hold on
plot(t_num_full, theta_res_full(2,:))
plot(t_num_full, lambda_res_full(2,:))
for ii=1:length(t_num)
    xline(t_num(ii),'k:')
end
hold off
subplot(4,1,3)
hold on
plot(t_num_full, theta_res_full(3,:))
plot(t_num_full, lambda_res_full(3,:))
for ii=1:length(t_num)
    xline(t_num(ii),'k:')
end
hold off
subplot(4,1,4)
hold on
plot(t_num_full, theta_res_full(4,:))
plot(t_num_full, lambda_res_full(4,:))
for ii=1:length(t_num)
    xline(t_num(ii),'k:')
end
hold off


figure
subplot(2,1,1)
hold on
plot(t_num_full, x_res_full(1,:))
for ii=1:length(t_num)
    xline(t_num(ii),'k:')
end
hold off
subplot(2,1,2)
hold on
plot(t_num_full, x_res_full(2,:))
for ii=1:length(t_num)
    xline(t_num(ii),'k:')
end
hold off
