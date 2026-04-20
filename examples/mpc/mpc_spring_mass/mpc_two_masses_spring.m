%% Minimal MPC example for algorithm testing with many switching events.
% This script compares real-time MPC and fully converged MPC variants for a
% piecewise-smooth two-mass spring system and can optionally simulate
% plant-model mismatch using a separate integrator.
%%
clear all; clc; close all;

import casadi.*
import nosnoc.*

use_rtmpc = 1; % Use real-time MPC algorithms instead of fully converged MPC with warm starts
use_ccopt = 0; % Recommended: Use our fastest QPCC/MPCC solver; requires installation (see nosnoc README).
save_qpecs = 0; % Save QPCC subproblems
model_plant_mismatch = 1; % MPC and simulation model use different integrators/models
visualize = 1; 
%% 
N_stages = 10; % Number of control intervals in the MPC problem
T = 1.0; % MPC prediction horizon
N_steps = 60; % Number of MPC simulation steps
use_fesd = 1;  % Use the FESD discretization method
N_FE = 2; % Number of integration steps per control interval

%% model
% Define variables
p1 = SX.sym('p1');  % Position of mass 1
p2 = SX.sym('p2');  % Position of mass 2
v1 = SX.sym('v1');  % Velocity of mass 1
v2 = SX.sym('v2');  % Velocity of mass 2
u  = SX.sym('u');   % Control input

x = [p1; p2; v1; v2];  % State vector

% Parameters
m1 = 1.0;  % Mass 1
m2 = 1.0;  % Mass 2
k  = 15.0; % Spring stiffness
Fc = 1.0;  % Coulomb friction force

% ODE without friction
dp1 = v1;
dp2 = v2;
dv1 = (1/m1) * (u - k * (p1 - p2));
dv2 = (1/m2) * (k * (p1 - p2));

f_base =  [dp1; dp2; dv1; dv2];
friction_1 = [0;0;- (1/m1) * Fc * sign(v1);0];
friction_2 = [0;0; 0;- (1/m2) * Fc * sign(v2)];

x0 = [-1.0; +1.0; 0.5; -0.5];  % [p1, p2, v1, v2]

Q = diag([10, 10, 1, 1]); % Penalize position errors more than velocity errors
R = 1.0;  % Control effort penalty
x_ref = [0; 0; 0; 0];  % Reference state
u_ref = 0;

% Objective function
f_q = (x - x_ref)' * Q * (x - x_ref) + (u - u_ref)' * R * (u - u_ref);
f_q_T = 10*(x - x_ref)' * Q * (x - x_ref);

u_ref = 0;  % Desired control input
x_ref = [0; 0; 0; 0];  % Reference state
% Cost function weights
Q = diag([10, 10, 1, 1]); % Penalize position errors more than velocity errors
R = 0.1;  % Control effort penalty

% Bounds
lbx = -inf(4,1); % No lower bounds on states
ubx =  inf(4,1); % No upper bounds on states
lbu = -5.0; % Lower bound on control
ubu =  5.0; % Upper bound on control

%% all options
problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
homotopy_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosnoc problems.
mpecopt_options = mpecopt.Options(); % Initialize all options related to the MPEC solver used for solving nosnoc problems.
ccopt_options = nosnoc.ccopt.Options(); % CCopt options
mpc_options = nosnoc.mpc.Options(); % Initialize all options related to the MPC implementation.

ccopt_options.opts_madnlp.linear_solver = 'Ma27Solver';
ccopt_options.opts_ccopt.relaxation_update.TYPE = 'RolloffRelaxationUpdate';
ccopt_options.opts_ccopt.relaxation_update.rolloff_slope = 2.0;
ccopt_options.opts_ccopt.relaxation_update.rolloff_point = 1e-6;
ccopt_options.opts_ccopt.relaxation_update.sigma_min = 1e-6;
ccopt_options.opts_madnlp.tol=1e-6;
ccopt_options.opts_ccopt.q_regularization = 'critical_rho';
ccopt_options.opts_ccopt.critical_rho_factor = 0.9999;

% Choose the Runge-Kutta method and number of stages
problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
problem_options.n_s = 2; % Number of stage points in the RK method (determines accuracy)
problem_options.cross_comp_mode = "FE_FE";
homotopy_options.homotopy_steering_strategy = "DIRECT";

% Time settings
problem_options.N_stages = N_stages; % Number of control intervals over the time horizon.
problem_options.N_finite_elements = N_FE; % Number of finite elements (integration steps) per control interval (optionally a vector may be passed).
problem_options.T = T;    % Time horizon (for a time-optimal problem, the actual horizon length is obtained by solving the problem).
problem_options.use_fesd = use_fesd;
%% nosnoc models
model = nosnoc.model.Pss(); %
model.x = x;
model.x0 = x0;
model.lbx = lbx; % Lower bounds on states
model.ubx = ubx; % Upper bounds on states
% define control vectors
model.u = u;
model.lbu = lbu;
model.ubu = ubu;

f_11 = [0;0;+(1/m1)*Fc;0];
f_12 = [0;0;-(1/m1)*Fc;0];

f_21 = [0;0;0;+(1/m2)*Fc];
f_22 = [0;0;0;-(1/m2)*Fc];

c1 = v1;
c2 = v2;
% Sign matrix for the modes
S1 = [-1;1];
S2 = [-1;1];
F1 = [f_11 f_12];
F2 = [f_21 f_22];

model.S = {S1,S2};
model.c = {c1,c2};
model.F = {F1 F2};
model.f_0 = f_base;

model.f_q = f_q;
model.f_q_T = f_q_T;
% 
% model.g_path = p1-p2;
% model.lbg_path = -inf;
% model.ubg_path = 0;

%% Setup solver an mpc options
homotopy_options.complementarity_tol = 1e-6;
homotopy_options.N_homotopy = 7;
homotopy_options.lift_complementarities = 0;
%solver_options.print_level = 0;
homotopy_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; % Using an HSL solver is significantly faster than the default 'mumps', but requires installation
% homotopy_options.opts_casadi_nlp.ipopt.hsllib = '/home/anton/tools/HSL_jll.jl-2023.11.7/override/lib/x86_64-linux-gnu-libgfortran5/libhsl.so';

%% QPCC/QPEC solver options
qpcc_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosnoc problems.
% qpcc_options.homotopy_update_rule = 'superlinear';
qpcc_options.homotopy_update_slope = 0.1;
qpcc_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; % Using an HSL solver is significantly faster than the default 'mumps', but requires installation
% qpcc_options.opts_casadi_nlp.ipopt.hsllib = '/home/anton/tools/HSL_jll.jl-2023.11.7/override/lib/x86_64-linux-gnu-libgfortran5/libhsl.so';
qpcc_options.homotopy_steering_strategy = "DIRECT";
qpcc_options.opts_casadi_nlp.ipopt.tol = 1e-9;
qpcc_options.complementarity_tol  = 1e-6;
qpcc_options.print_level  = 5;

gurobi_options = nosnoc.qpec.GurobiOptions();
gurobi_options.method = 'reg';
% gurobi_options.method = 'miqp';
% gurobi_options.lower_bounds_comps = false;
gurobi_options.method = 'sos1';

%% Mpc
mpc_options.do_shift_initialization = false;
mpc_options.warmstart_qpec = true;
mpc_options.solve_advanced_problem = true;
mpc_options.advanced_n_qpecs = 1;
mpc_options.fast_sigma_0 = 1e0;
mpc_options.discard_constraints_in_hessian = false;
mpc_options.sqpec_hessian_convexification = "MIRROR";
mpc_options.on_qpec_failure = "keep";
%% create mpc object
if use_rtmpc
    if use_ccopt
        mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, ccopt_options, ccopt_options);
    else
        mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, homotopy_options, homotopy_options);
    end
else
    if use_ccopt
        mpc = nosnoc.mpc.FullMpcc(model, mpc_options, problem_options, ccopt_options);
    else
        mpc = nosnoc.mpc.FullMpcc(model, mpc_options, problem_options, homotopy_options);
    end
end

%% Create sim model and integrator

% Create the true model with air resistance
sim_model = nosnoc.model.Pss(); % Initialize a nosnoc model (depending on the underlying nonsmooth model, several options exist)
% define differential states and populate the model.
sim_model.x = x;
sim_model.x0 = x0;
sim_model.lbx = -inf(4,1);
sim_model.ubx = inf(4,1);
% define control vectors
sim_model.u = u;
sim_model.lbu = -inf(1,1);
sim_model.ubu = inf(1,1);
% Dynamics of the piecewise smooth systems


% Define the regions of the PSS
sim_model.c = {c1 c2};
sim_model.S = {S1 S2};
sim_model.F = {F1 F2}; % The columns of this matrix store the vector fields for each region.
sim_model.f_0 = f_base;

sim_problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
integrator_options = nosnoc.integrator.Options();
sim_solver_options = integrator_options.fesd_solver_opts; % Initialize all options related to the MPEC solver used for solving nosnoc problems.

% Choose the Runge-Kutta method and number of stages
sim_problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
sim_problem_options.n_s = 3; % Number of stage points in the RK method (determines accuracy)

% Time settings
sim_problem_options.N_finite_elements = 2; % Number of finite elements (integration steps) per control interval (optionally a vector may be passed).
sim_problem_options.T_sim = problem_options.h;
sim_problem_options.N_sim = 5;
sim_problem_options.print_level = 0;

% Simulation solver options
sim_solver_options.print_level = 0;
sim_solver_options.homotopy_steering_strategy = HomotopySteeringStrategy.ELL_INF; % Use the $\ell_{\infty}$ steering strategy
sim_solver_options.complementarity_tol = 1e-10; % Value used to drive the complementarity residual to zero.
% sim_solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; % Using an HSL solver is significantly faster than the default 'mumps', but requires installation
integrator = nosnoc.Integrator(sim_model, sim_problem_options, integrator_options);
%% Do MPC with accurate model

plot_intermediate_solutions = true;
if plot_intermediate_solutions
    q_plot = subplot(311); hold on;
    xlim([0 (N_steps+N_stages)*problem_options.h])
    ylim([-pi*1.1 pi*1.1])
    xlabel(q_plot, '$t$')
    ylabel(q_plot, '$p$')

    v_plot = subplot(312); hold on;
    xlim([0 (N_steps+N_stages)*problem_options.h])
    ylim([-5 5])
    xlabel(v_plot, '$t$')
    ylabel(v_plot, '$v$')

    u_plot = subplot(313); hold on;
    xlim([0 (N_steps+N_stages)*problem_options.h])
    ylim([-11 11])
    xlabel(u_plot, '$t$')
    ylabel(u_plot, '$u$')
end

x = model.x0; u = []; t = 0; tf = [];
x0 = x;
last_step = N_steps;
for step=1:last_step%N_steps
    [u_i, stats] = mpc.get_feedback(x0);
    tf_i = stats.feedback_time;
    tf = [tf, tf_i];
    fprintf("MPC step: %d, Feedback time: %d\n", step, tf_i);
    if step ~= last_step
        mpc.do_preparation();
    else
        break
    end
    if model_plant_mismatch
        integrator.set_x0(x0);
        [t_grid, x_sim] = integrator.simulate("u", repmat(u_i, [1,5]), "x0", x0);
        x0 = x_sim(:, end);
    else
        x0 = mpc.get_predicted_state();
    end
    if use_rtmpc
        if save_qpecs && mod(step-1, 5) == 0
            save_qpec(mpc.qpec, ['SPRING_FRICTION_CONVEX_', num2str(step)])
        end
    end

    % Plot intermediate solution by getting x, t, and u from mpc object
    % Using `.get`
    if plot_intermediate_solutions
        x_res = mpc.get('x');
        q_res = x_res(1:2,:);
        v_res = x_res(3:4,:);
        t_grid = mpc.get_time_grid();
        u_res = mpc.get('u');
        t_grid_u = mpc.get_control_grid();

        %
        cla(q_plot); hold on;
        cla(v_plot); hold on;
        cla(u_plot); hold on;

        plot(q_plot, t, x(1:2,:))
        plot(q_plot, t(end)+t_grid, q_res)
        yline(q_plot,x_ref(1:2),'k--')
        xlabel(q_plot, '$t$')
        ylabel(q_plot, '$p$')

        plot(v_plot, t, x(3:4,:))
        plot(v_plot, t(end)+t_grid, v_res)
        yline(v_plot,x_ref(3:4),'k--')
        xlabel(v_plot, '$t$')
        ylabel(v_plot, '$v$')

        if step>1
            stairs(u_plot, t, [u, u(:,end)]')
        end
        stairs(u_plot, t_grid_u+t(end), [u_res, u_res(:,end)]')
        yline(u_plot,u_ref,'k--')
        xlabel(u_plot, '$t$')
        ylabel(u_plot, '$u$')
% 
    end

    x = [x, x0];
    u = [u,u_i];
    t = [t, t(end) + problem_options.h];
end

%% Plot
figure
latexify_plot()
subplot(311)
plot(t,x(1:2,:),'LineWidth',1.5)
hold on;
yline(x_ref(1:2),'k--','LineWidth',1.5)
xlabel("$t$")
ylabel("$q$")
subplot(312)
plot(t,x(3:4,:),'LineWidth',1.5)
hold on;
yline(x_ref(3:4),'k--','LineWidth',1.5)
xlabel("$t$")
ylabel("$v$")
% ylim([-5 25])
subplot(313)
stairs(t,[u,u(:,end)]','LineWidth',1.5)
hold on
% yline(u_max,'k--','LineWidth',1.5)
% yline(-u_max,'k--','LineWidth',1.5)
xlabel("$t$")
ylabel("$u$")
% ylim([-1.1*u_max 1.1*u_max])


%% Visualize

% Given state history X (from MPC simulation)
if visualize
figure;
hold on;
grid on;
xlabel('Position');
ylabel('Masses');
title('Oscillation of Two Masses with Spring');

% Define fixed axis limits
x_min = min(x(:,1:2), [], 'all') - 0.5; % Slight margin
x_max = max(x(:,1:2), [], 'all') + 0.5;
num_points = 20;  % Number of wave segments
% Animation loop
num_steps = size(x, 2);
for i = 1:num_steps
    clf;
    hold on;
    grid on;
    
    % Mass positions
    p1 = x(1,i);
    p2 = x(2,i);
    
    % Plot masses
    plot(p1, 0, 'ro', 'MarkerSize', 15, 'MarkerFaceColor', 'r'); % Mass 1
    plot(p2, 0, 'bo', 'MarkerSize', 15, 'MarkerFaceColor', 'b'); % Mass 2
    
    % Plot wavy spring (discretized sinusoidal curve)
    
    spring_x = linspace(p1, p2, num_points);
    spring_y = 0.05 * sin(linspace(0, 4*pi, num_points)); % Wave shape
    
    plot(spring_x, spring_y, 'k-', 'LineWidth', 2); % Draw spring

    % Set fixed axis limits
    xlim([x_min, x_max]);
    ylim([-0.2, 0.2]);  
    axis equal;
    
    drawnow;
    pause(0.25);  % Adjust animation speed
end
end