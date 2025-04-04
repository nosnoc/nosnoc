clear;clc;close all;
import casadi.*
import nosnoc.*


%% populate options
problem_options = nosnoc.Options(); % problem_options = nosnoc.Options();
solver_options = nosnoc.reg_homotopy.Options();

%% set some options
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.rk_representation = RKRepresentation.integral;
problem_options.n_s = 2;
problem_options.N_stages = 20; % number of control intervals
problem_options.N_finite_elements = 2; % number of finite element on every control interval (optionally a vector might be passed)
problem_options.T = 15;    % Time horizon
problem_options.dcs_mode = "Stewart"; % or "Heaviside"
problem_options.use_numerical_clock_state = true;

problem_options.step_equilibration = StepEquilibrationMode.heuristic_mean;

%problem_options.x_box_at_stg = 0;
%problem_options.x_box_at_fe = 0;
problem_options.relax_terminal_constraint = ConstraintRelaxationMode.ELL_2;
%problem_options.relax_terminal_numerical_time = false;

solver_options.homotopy_steering_strategy = "ELL_1";
solver_options.decreasing_s_elastic_upper_bound = true;
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
solver_options.print_level=5;
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
% lsq on control
model.lsq_u = {u, 0, 1};

%% create solver and solve problem
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.set_param('rho_terminal', 100);
ocp_solver.solve();
%% get results and plot
x = ocp_solver.get('x');
u = ocp_solver.get('u');
%theta = ocp_solver.get('theta');
t_grid = ocp_solver.get_time_grid();
t_grid_u = ocp_solver.get_control_grid();

figure
plot(t_grid, x);
figure
stairs(t_grid_u, [u,u(end)])
% figure
% plot(t_grid(2:end), theta);
