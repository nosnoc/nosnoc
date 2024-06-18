clear;clc;close all;
import casadi.*
import nosnoc.*


%% populate options
problem_options = nosnoc.Options(); % problem_options = NosnocProblemOptions();
solver_options = nosnoc.solver.Options();
%% set some options
problem_options.rk_scheme = RKSchemes.GAUSS_LEGENDRE;
problem_options.rk_representation = RKRepresentation.differential_lift_x;
problem_options.n_s = 2;
problem_options.N_stages = 20; % number of control intervals
problem_options.N_finite_elements = 2*ones(20,1); % number of finite element on every control interval (optionally a vector might be passed)
problem_options.N_finite_elements(1:5) = 4;
% problem_options.N_finite_elements(10:end) = 1;
problem_options.T = 20;    % Time horizon
problem.options.dcs_mode = "Stewart"; % or "Heaviside"

problem_options.x_box_at_stg = 0;
problem_options.x_box_at_fe = 0;
problem_options.relax_terminal_constraint = "NONE";
problem_options.rho_terminal = 1e4; 
%problem_options.relax_terminal_numerical_time = false;
%% Create model
model = nosnoc.model.Pss(); 
q = SX.sym('q'); 
v = SX.sym('v'); 
u = SX.sym('u');
model.x = [q;v];
model.x0 = [0;0];
model.lbx = [-inf;-20]; model.ubx = [inf;20];
model.u = u;
model.lbu = -10; model.ubu = 10;
% Dyanmics and the regions
f_1 = [v;u]; % mode 1 - nominal
f_2 = [v;3*u]; % mode 2 - turbo
model.S = [-1;1];
model.F = [f_1 f_2];
model.c = v-10;
% Add terminal constraint
model.g_terminal = [q-200;v-0];
% lsq on control
model.lsq_u = {u, 0, 1};

%% create solver and solve problem
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();

x = ocp_solver.get('x');
u = ocp_solver.get('u');
t_grid = ocp_solver.get_time_grid();
t_grid_u = ocp_solver.get_control_grid();
%% plot
figure
subplot(121)
plot(t_grid, x(1,:));
grid on
subplot(122)
plot(t_grid, x(2,:));
grid on
figure
stairs(t_grid_u, [u,u(end)])
% how to create an integrator?
% integrator = nosnoc.Integrator(model, problem_options, solver_options, [], []); % What could be further optional argumetns, i would prefer a varargin instead of passing empty stuff.
% [results,stats] = integrator.solve();
