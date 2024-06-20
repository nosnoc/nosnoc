clear;clc;close all;
import casadi.*
import nosnoc.*


%% populate options
problem_options = nosnoc.Options(); % problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();

%% set some options
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 2;
problem_options.time_optimal_problem = 1;
problem_options.N_stages = 11; % number of control intervals
problem_options.N_finite_elements = 2; % number of finite element on every control interval (optionally a vector might be passed)
problem_options.T = 1;    % Time horizon

problem_options.step_equilibration = StepEquilibrationMode.linear_complementarity;

% problem.options.dcs_mode = "Stewart"; % or "Heaviside"
%% Create model
model = nosnoc.model.Heaviside();  % same as PSS+Heaviside, this creates a nosnoc.dcs.heaviside
x = SX.sym('x');
u = SX.sym('u');
alpha = SX.sym('alpha'); % step function 
c = x-1;
model.x = x;
model.x0 = 0;
model.lbx = [-inf]; model.ubx = [inf];
model.u = u;
model.lbu = -5; model.ubu = 5;
model.alpha = alpha;
% Dyanmics and the regions
f_1 = -1+0.05*x+u;
f_2 = 1+0.05*x+u;
f_rhs = alpha*f_1 + (1-alpha)*f_2; % more complicated expressions and with alpha being a vector are possible as well;
model.f_x = f_rhs;
model.c = c;
% Add terminal constraint
model.g_terminal = [x-2];

%% create solver and solve problem
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();

x = ocp_solver.get('x');
u = ocp_solver.get('u');
t_grid = ocp_solver.get_time_grid();
t_grid_u = ocp_solver.get_control_grid();

figure
plot(t_grid, x);
figure
stairs(t_grid_u, [u,u(end)])
