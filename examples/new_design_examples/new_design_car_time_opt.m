clear;clc;close all;
import casadi.*
import nosnoc.*


%% populate options
problem_options = nosnoc.Options(); % problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();

%% set some options
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.rk_representation = RKRepresentation.integral;
problem_options.n_s = 2;
problem_options.N_stages = 11; % number of control intervals
problem_options.N_finite_elements = 4; % number of finite element on every control interval (optionally a vector might be passed)
problem_options.T = 1.0;    % Time horizon
problem_options.dcs_mode = "Heaviside"; % or "Heaviside"
problem_options.time_optimal_problem = true;
problem_options.step_equilibration = StepEquilibrationMode.linear_complementarity;
problem_options.use_fesd = 1;
problem_options.use_speed_of_time_variables = 1;
problem_options.local_speed_of_time_variable = 1;
%% Create model
model = nosnoc.model.Pss(); 
q = SX.sym('q'); 
v = SX.sym('v'); 
u = SX.sym('u');
model.x = [q;v];
model.x0 = [0;0];
model.lbx = [-inf;-20]; model.ubx = [inf;20];
model.u = u;
model.u0 = 1;
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

%% Plot resukts
x = ocp_solver.get('x');
u = ocp_solver.get('u');
t_grid = ocp_solver.get_time_grid();
t_grid_u = ocp_solver.get_control_grid();
figure
plot(t_grid, x(1,:));
figure
plot(t_grid, x(2,:));
figure
stairs(t_grid_u, [u,u(end)])

if problem_options.use_fesd 
    h_res = ocp_solver.discrete_time_problem.w.h.res;
    figure
    stairs(t_grid, [h_res, h_res(end)])
end
% how to create an integrator?
% integrator = nosnoc.integrator.FESD(model, problem_options, solver_options, [], []); % What could be further optional argumetns, i would prefer a varargin instead of passing empty stuff.
% [t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();
