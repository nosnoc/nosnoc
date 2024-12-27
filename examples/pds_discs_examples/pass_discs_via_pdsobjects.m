% Optimal control example from 
% Finite Elements with Switch Detection for Numerical Optimal Control of
% Projected Dynamical Systems, Proceedings of the IEEE Conference on Decision and Control (CDC) 2024
% Anton Pozharskiy, Armin Nurkanovic, Moritz Diehl

% *** This is the same example as ''pass_discs.m`` but now useding the PDSObjects
% class for a more user friendly and modular problem defintion ****
%%
clear all
close all
import casadi.*

%% create nosnoc model and options objects
% model = nosnoc.model.Pds();
model = nosnoc.model.PDSObjects; % Initialize model which is a container for the objects.
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();

%% Parameters
T = 5;  % time horizon
R = 1; % radius of manipulators
R_obj = 2; % radius of manipulated object
n_d = 2; % number of dimensions;
%% Define projected system
% control "forces" for manipulators
u1 = SX.sym('u1', n_d);
u2 = SX.sym('u2', n_d);
u = [u1;u2];

model.u = u;
model.u0 = [0;0;0;0];
model.lbu = [-10;-10;-10;-10];
model.ubu = [10;10;10;10];

% Define objects
ball1 = nosnoc.objects.Ball(R,n_d); % Manipulator 1
ball2 = nosnoc.objects.Ball(R,n_d); % Manipulator 2
ball3 = nosnoc.objects.Ball(R_obj,n_d); % Manipulated object 


% Define their initial state, bounds and dynamics
ball1.x0 = [-10;10]; % Initial position of ball
ball1.x_dot = u1; % Ball is actuated directly with the controls.
ball1.lbx = [-inf;-inf];
ball1.ubx = [-2.5;inf];
ball1_target = [-10;0];

ball2.x0 = [5;-3]; % Initial position of ball
ball2.x_dot = u2; % Ball is actuated directly with the controls.
ball2.lbx = [2.5;-inf];
ball2.ubx = [inf;inf];
ball2_target = [10;0];

ball3.x0 = [0;0]; % Initial position of ball
ball3.x_dot = [0;0]; % Ball is not directly actuated
ball3.lbx = [-inf;-inf];
ball3.ubx = [inf;inf];
ball3_target = [0;10];


x_target = vertcat(ball1_target, ball2_target, ball3_target); % Total target state.
% Contacts 
model.addContact(ball1, ball3); % Add contact between manipulaor 1 and obj
model.addContact(ball2, ball3); % Add contact between manipulaor 2 and obj
% costs
model.f_q = 1e-4*norm_2(model.u)^2;
model.f_q_T = (model.x-x_target)'*diag([1e0,1e0,1e0,1e0,1e3,1e3])*(model.x-x_target);

% Time discertization settings
problem_options.T = T;
problem_options.N_stages = 30; % numbe of control stages\intervals
problem_options.N_finite_elements = 3; % number of integration steps in every control interval
problem_options.n_s = 2;

% Homotopy solver options
solver_options.complementarity_tol = 1e-10;
solver_options.print_level = 3;
solver_options.opts_casadi_nlp.ipopt.max_iter = 500;
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; % if solver fails (i.e. nothing happens) because
% HLS library not instaled use 'mumps' instead of ma27

%% create solver and solve optimal control problem
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();
%% plot
x_res = ocp_solver.get("x");
u_res = ocp_solver.get("u");
h_res = ocp_solver.get("h");
t_res = ocp_solver.get_time_grid();
fig = figure('Position', [10 10 1600 800]);
rectangle('Position',[-1.5 -25 3 50],'FaceColor',[1 0 0 0.2],'LineStyle', '--', 'EdgeColor' , 'red')
rectangle('Position',[-2 8 4 4], 'Curvature', 1, 'FaceColor',[0.8500 0.3250 0.0980,0.5], 'LineStyle', '--', 'EdgeColor' , [0.8500 0.3250 0.0980])
rectangle('Position',[-11 -1 2 2], 'Curvature', 1, 'FaceColor',[0 0.4470 0.7410,0.5], 'LineStyle', '--', 'EdgeColor' , [0 0.4470 0.7410])
rectangle('Position',[9 -1 2 2], 'Curvature', 1, 'FaceColor',[0 0.4470 0.7410,0.5], 'LineStyle', '--', 'EdgeColor' , [0 0.4470 0.7410])
plot_pass_discs(h_res,x_res,[R,R,R_obj], ["circle","circle","circle"], fig, 'coop');
