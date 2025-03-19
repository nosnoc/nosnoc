% Optimal control example from 
% Finite Elements with Switch Detection for Numerical Optimal Control of
% Projected Dynamical Systems, Proceedings of the IEEE Conference on Decision and Control (CDC) 2024
% Anton Pozharskiy, Armin Nurkanovic, Moritz Diehl
%%
clear all
close all
import casadi.*

%% create nosnoc model and options objects
model = nosnoc.model.Pds();
problem_options = nosnoc.Options();
solver_options = nosnoc.reg_homotopy.Options();

%% parameter
T = 5;
R = 1; % radius of manipulators
R_obj = 1; % radius of manipulated object
R_obstacle = 5; % radius of obstacle
%% Define projected system
x1 = SX.sym('x1', 2); % center of manipulator 1
x2 = SX.sym('x2', 2); % center of manipulator 2
x3 = SX.sym('x3', 2); % center of manipulated object
gate = SX.sym('gate', 1); % vertical position of gate, horizontal fix
x = [x1;x2;x3;gate];
x_target = [-10;3;-10;0;-7;3;0];
x0 = [-10;2;10;4;7;3.5;0.9];
u1 = SX.sym('u1', 2);
u2 = SX.sym('u2', 2);

% populate nosnoc PDS model
model.x = x;
model.lbx = [-inf;-inf;-inf;-inf;-inf;-inf;-inf];
model.ubx = [inf;5;inf;5;inf;5;inf];
model.x0 = x0;
model.u = [u1;u2];
model.lbu = [-10;-10;-10;-10];
model.ubu = [10;10;10;10];
model.u0 = [0;0;0;0];
model.c = [norm_2(x3-x1)-(R+R_obj);
    norm_2(x3-x2)-(R+R_obj);
    norm_2(x2-x1)-(R+R);
    x1(2)-gate-R]; % these are the distance functions between manipualtors and objects
model.g_path = [x2(2)-gate-R;
    x3(2)-gate-R_obj;
    norm_2(x1-[0;5]) - R-R_obstacle;
    norm_2(x2-[0;5]) - R-R_obstacle;
    norm_2(x3-[0;5]) - R_obj-R_obstacle]; % path constraints for obstacle avoidance and not touching the gate
model.lbg_path = [0;0;0;0;0];
model.ubg_path = [inf;inf;inf;inf;inf];
model.f_x_unconstrained = [u1;u2;0;0;0]; % dynamics of the PDS when no constraints are active

% costs
model.f_q = 1e-4*norm_2(model.u)^2;
model.f_q_T = (x-x_target)'*diag([1,1,1,1,1e3,1e3,0])*(x-x_target);
%model.g_T = x3 - [-7;0];

% Time discertization settings

problem_options.T = T;
problem_options.N_stages = 30;
problem_options.N_finite_elements = 3;
problem_options.n_s = 2;

% Homotopy solver settings 
solver_options.complementarity_tol = 1e-10;
solver_options.print_level = 3;
solver_options.opts_casadi_nlp.ipopt.max_iter = 500;
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';

%% create solver
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();
%% plot
x_res = ocp_solver.get("x");
u_res = ocp_solver.get("u");
h_res = ocp_solver.get("h");
t_res = ocp_solver.get_time_grid();
fig = figure('Position', [10 10 1600 800]);
rectangle('Position',[-5, 0 10 10], 'Curvature', 1, 'FaceColor',[1 0 0 0.1], 'LineStyle', '--', 'EdgeColor' , [1 0 0])
rectangle('Position',[-8 2 2 2], 'Curvature', 1, 'FaceColor',[0.8500 0.3250 0.0980,0.5], 'LineStyle', '--', 'EdgeColor' , [0.8500 0.3250 0.0980])
rectangle('Position',[-11 -1 2 2], 'Curvature', 1, 'FaceColor',[0 0.4470 0.7410,0.5], 'LineStyle', '--', 'EdgeColor' , [0 0.4470 0.7410])
rectangle('Position',[-11 2 2 2], 'Curvature', 1, 'FaceColor',[0 0.4470 0.7410,0.5], 'LineStyle', '--', 'EdgeColor' , [0 0.4470 0.7410])
plot_open_gate(h_res,x_res,[R,R,R_obj], ["circle","circle","circle"], fig, 'open_gate');
