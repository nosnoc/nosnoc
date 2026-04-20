% This script gives a minimal introduction to MPC in nosnoc by setting up a 
% simple piecewise smooth vehicle model, defining an OCP with stage and 
% terminal costs, and running the controller in closed loop over several MPC 
% steps. It shows the main workflow: create the model, choose solver and MPC 
% options, compute feedback, update parameters, and inspect the predicted 
% trajectories.
% 
% For more advanced setups, it is worth looking at the other nosnoc examples, 
% which cover additional modeling and solver features. 
% For faster MPC execution, also consider the CCOpt-based solver variants.

%% Init
clear; clc; close all;
import casadi.*
%% Basic options
% Load default options for various nosnoc classes 
problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
reg_solver_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosonc problems.
mpc_options = nosnoc.mpc.Options(); % Initialize all options related to the MPC implementation. 
use_hyrti = true; % If true, use the Hybrid Real-Time Iteration, if false, solve OCP to convergence.

% Choosing the Runge - Kutta Method and number of stages
problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
problem_options.n_s = 2; % Number of stage points in the RK method (determines accuracy)
problem_options.cross_comp_mode = 'FE_FE';

% Time-settings  - Solve a 
problem_options.N_stages = 5; % Number of control intervals over the time horizon.
problem_options.N_finite_elements = 2; % Number of finite element (integration steps) on every control interval (optionally a vector might be passed).
problem_options.T = 5;    % Time horizon (for a time-optimal problem, the actual horizon length is obtained by solving the problem).

%% Create a nosnoc model
model = nosnoc.model.Pss(); % Initialize a nosnoc model, (depending on the underlying nonsmooth model, several options exist)
% bound parameters
v_max = 30;
u_max = 5;
% define differential states and populate the model.
q = SX.sym('q'); 
v = SX.sym('v'); % CasADi symbolic variables for states

model.x = [q;v]; % populate model state vectors
model.x0 = [0;0]; % initial value
model.lbx = [-inf;-v_max]; % lower bounds on states
model.ubx = [inf;v_max]; % upper bounds on states

v_target = SX.sym('v_target'); % CasADi symbolic variables for target v in running cost.
model.p_global = v_target;
model.p_global_val = 15;

% define control vectors
u = SX.sym('u');  % CasADi symbolic variables for controls
model.u = u;
model.lbu = -u_max; 
model.ubu = u_max;

%% Dynamics of the piecewise smooth systems
f_1 = [v;u]; % mode 1 - nominal
f_2 = [v;3*u]; % mode 2 - turbo
% Define the regions of the PSS
v_threshold = 10;
model.c = v-v_threshold; % single switching functions (implies two regions)
model.S = [-1;1]; % Region R_1 is defined by c<0 (hence the -1), region R_2 by c>0 (hence the +1) in the sign matrix S.
model.F = [f_1 f_2]; % The columns of this matrix store the vector fields of every region.

q_goal = 400;
v_goal = 0;
model.f_q = (q-q_goal)^2 + u^2 +(v-v_target)^2; % Add stage cost
model.f_q_T = 10*((q-q_goal)^2 + (v-v_goal)^2); % Add terminal quadratic cost

N_steps = 30; % number of MPC steps;

% Setup solver an mpc options
reg_solver_options.homotopy_update_rule = 'superlinear'; % Use superlinear update rule for relaxation parameter sigma.
reg_solver_options.homotopy_update_slope = 0.05; % Rate which the relaxation sigma is reduced at: sigma_i+1 = kappa*sigma_i 
reg_solver_options.homotopy_update_exponent = 2; % Rate which the relaxation sigma is reduced at: sigma_i+1 = sigma_i^kappa
reg_solver_options.complementarity_tol = 1e-7; % Value to drive the complementarity residual to.
reg_solver_options.N_homotopy = 10; % Maximum number of homotopy iterations.
reg_solver_options.print_level = 0;
% homotopy_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; % Using an HSL solver is a significant speedup over the default 'mumps' but requires installation
% Set the fast sigma for followup solves to the complementarity tol to only do 1 nlp solve. 
mpc_options.fast_sigma_0 = 1e-7;

% create mpc object
if use_hyrti
    reg_solver_options.homotopy_update_rule = 'linear'; % Use superlinear update rule for relaxation parameter sigma.
    reg_solver_options.homotopy_update_slope = 0.1; % Rate which the relaxation sigma is reduced at: sigma_i+1 = kappa*sigma_i 
    mpc_options.fast_sigma_0 = 1e0;
    mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, reg_solver_options, reg_solver_options);
else
    mpc = nosnoc.mpc.FullMpcc(model, mpc_options, problem_options, reg_solver_options);
end

% Flag for plotting intermediate solutions
plot_intermediate_solutions = true;
if plot_intermediate_solutions
    q_plot = subplot(311);
    v_plot = subplot(312);
    u_plot = subplot(313);
end

% Do MPC assuming the predicted state is accurate, in practice this may not be true.
x = model.x0; u = []; t = 0; tf = [];
x0 = x;

% Calculate target v
v_target_val = (x0(1) - q_goal)/problem_options.T;
mpc.set_param('p_global', {}, v_target_val) % for p_global the index parameter should be empty, otherwise it should be the control stage for p_time_var.
for step=1:N_steps
    [u_i, stats] = mpc.get_feedback(x0);
    tf_i = stats.feedback_time;
    tf = [tf, tf_i];
    fprintf("MPC step: %d, Feedback time: %d\n", step, tf_i);
    mpc.do_preparation();
    x0 = mpc.get_predicted_state();

    % Plot intermediate solution by getting x, t, and u from mpc object
    % Using `.get`
    if plot_intermediate_solutions
        x_res = mpc.get('x');
        q_res = x_res(1,:);
        v_res = x_res(2,:);
        t_grid = mpc.get_time_grid();
        u_res = mpc.get('u');
        t_grid_u = mpc.get_control_grid();

        plot(q_plot, t_grid, q_res)
        plot(v_plot, t_grid, v_res)
        stairs(u_plot, t_grid_u, [u_res, u_res(end)])
    end 

    % Calculate target v
    v_target_val = (x0(1) - q_goal)/problem_options.T;
    % Set target velacity global parmeter.
    mpc.set_param('p_global', {}, v_target_val) % for p_global the index parameter should be empty, otherwise it should be the control stage for p_time_var.
    
    % x0_distrubance = (-2+4*rand(size(x0)));
    % x0 =  x0+x0_distrubance;
    x = [x, x0];
    u = [u,u_i];
    t = [t, t(end) + problem_options.h];
end

%% Plot
figure
latexify_plot()
subplot(311)
plot(t,x(1,:))
xlabel("$t$")
ylabel("$q$")
ylim([-5 q_goal+5])
subplot(312)
plot(t,x(2,:))
xlabel("$t$")
ylabel("$v$")
ylim([-5 25])
subplot(313)
stairs(t,[u,u(end)])
xlabel("$t$")
ylabel("$u$")
ylim([-1.1*u_max 1.1*u_max])

figure
stairs(tf)
ylabel("feedback time (s)")
xlabel("step")
