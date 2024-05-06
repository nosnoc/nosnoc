clear;clc;close all;
import casadi.*
import nosnoc.*


%% populate options
problem_options = nosnoc.core.Options(); % problem_options = NosnocProblemOptions();
problem_options = nosnoc.Options(); % problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();

%% set some options
problem_options.irk_scheme = IRKSchemes.RADAU_IIA;
problem_options.n_s = 2;
problem_options.time_optimal_problem = 1;
problem_options.N_stages = 10; % number of control intervals
problem_options.N_finite_elements = 2; % number of finite element on every control interval (optionally a vector might be passed)
problem_options.T = 1;    % Time horizon
% problem.options.dcs_mode = "Stewart"; % or "Heaviside"
%% Create model
% model = nosnoc.model.stewart();
model = nosnoc.model.heaviside();  % same as PSS+Heaviside, this creates a nosnoc.dcs.heaviside
x = SX.sym('x');
u = SX.sym('u');
alpha = SX.sym('alpha'); % step function 
c = x-1;
model.x = x;
model.x0 = 1;
model.lbx = [-inf]; model.ubx = [inf;20];
model.u = u;
model.lbu = -5; model.ubu = 5;
% Dyanmics and the regions
f_1 = -1+0.05*x+u;
f_2 = 1+0.05*x+u;
f_rhs = alpha*f_1 + (1-alpha)*f2; % more complicated expressions and with alpha being a vector are possible as well;
model.f_x = f_rhs;
model.c = c;
% Add terminal constraint
model.g_terminal = [q-200;v-0];

%% create solver and solve problem
mpcc = nosnoc.mpcc(problem_options, model); % WHY call here the mpcc creation? this should be internal
% mpcc = NosnocMPCC(problem_options, model);
% solver = NosnocSolver(mpcc, solver_options);
% solver = nosnoc.solver(mpcc, solver_options);
solver = nosnoc.solver(model, solver_options); % PREFERED OPTION
[results,stats] = solver.solve();


