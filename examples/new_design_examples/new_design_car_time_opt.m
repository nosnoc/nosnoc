clear;clc;close all;
import casadi.*
import nosnoc.*


%% populate options
problem_options = nosnoc.Options(); % problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();

%% set some options
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.rk_representation = RKRepresentation.integral;
problem_options.n_s = 2;
problem_options.N_stages = 11; % number of control intervals
problem_options.N_finite_elements = 2; % number of finite element on every control interval (optionally a vector might be passed)
problem_options.T = 1.0;    % Time horizon
problem.options.dcs_mode = "Stewart"; % or "Heaviside"
problem_options.time_optimal_problem = true;
problem_options.step_equilibration = StepEquilibrationMode.heuristic_mean;
%% Create model
% model = nosnoc.model.stewart();
model = nosnoc.model.Pss(); 
q = SX.sym('q'); 
v = SX.sym('v'); 
u = SX.sym('u');
model.x = [q;v];
model.x0 = [0;0];
model.lbx = [-inf;-20]; model.ubx = [inf;20];
model.u = u;
model.lbu = -5; model.ubu = 5;
% Dyanmics and the regions
f_1 = [v;u]; % mode 1 - nominal
f_2 = [v;3*u]; % mode 2 - turbo
model.S = [-1;1];
model.F = [f_1 f_2];
model.c = v-10;
% Add terminal constraint
model.g_terminal = [q-200;v-0];

%% create solver and solve problem
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();

x = ocp_solver.get_x();
u = ocp_solver.getU();
t_grid = ocp_solver.get_time_grid();
t_grid_u = ocp_solver.get_control_grid();
h_res = ocp_solver.discrete_time_problem.w.h.res;

figure
plot(t_grid, x);
figure
stairs(t_grid_u, [u,u(end)])
figure
stairs(t_grid, [h_res, h_res(end)])
% how to create an integrator?
% integrator = nosnoc.Integrator(model, problem_options, solver_options, [], []); % What could be further optional argumetns, i would prefer a varargin instead of passing empty stuff.
% [results,stats] = integrator.solve();
