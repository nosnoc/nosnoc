% Optimal control example from 
% Finite Elements with Switch Detection for Numerical Optimal Control of
% Projected Dynamical Systems, Proceedings of the IEEE Conference on Decision and Control (CDC) 2024
% Anton Pozharskiy, Armin Nurkanovic, Moritz Diehl
%%
clear all
close all
import casadi.*
%%
generate_video = false;
%% create nosnoc model and options objects
model = nosnoc.model.Pds();
problem_options = nosnoc.Options();
solver_options = nosnoc.reg_homotopy.Options();

%% Parameters
T = 5;  % time horizon
R = 1; % radius of manipulators
R_obj = 2; % radius of manipulated object
%% Define projected system
x1 = SX.sym('x1', 2); % center of manipulator disc 1
x2 = SX.sym('x2', 2); % center of manipulator disc 2
x3 = SX.sym('x3', 2); % center of manipulated disc
x = [x1;x2;x3];
x0 = [-10;10;5;-3;0;0];
x_target = [-10;0;10;0;0;10];
% control "forces" for manipulators
u1 = SX.sym('u1', 2);
u2 = SX.sym('u2', 2);

% Populate nosnoc PDS model
model.x = x;
model.lbx = [-inf;-inf;2.5;-inf;-inf;-inf];
model.ubx = [-2.5;inf;inf;inf;inf;inf];
model.x0 = x0;
model.u = [u1;u2];
model.lbu = [-10;-10;-10;-10];
model.ubu = [10;10;10;10];
model.u0 = [0;0;0;0];
model.c = [norm_2(x3-x1)-(R+R_obj);norm_2(x3-x2)-(R+R_obj)];
model.f_x_unconstrained = [u1;u2;0;0];

% costs
model.f_q = 1e-4*norm_2(model.u)^2;
model.f_q_T = (x-x_target)'*diag([1e0,1e0,1e0,1e0,1e3,1e3])*(x-x_target);

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
if generate_video
    plot_pass_discs(h_res,x_res,[R,R,R_obj], ["circle","circle","circle"], fig, 'coop');
else
    plot_pass_discs(h_res,x_res,[R,R,R_obj], ["circle","circle","circle"], fig);
end
