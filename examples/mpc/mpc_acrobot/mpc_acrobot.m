%% MPC example for an acrobot with joint friction.
% The MPC problem is to swing the system from the hanging configuration to
% the upright reference while penalizing state tracking error and control
% effort. The model is piecewise smooth because Coulomb friction at both
% joints changes sign with the joint velocities.

clear all; clc; close all;

import casadi.*
import nosnoc.*
%% Basic options
use_rtmpc = 1;  % Use real-time MPC instead of solving the problem to full convergence
save_qpecs = 0; % Save QPCC subproblems (for benchmarking)
use_ccopt = 0; % Recommended: Use our fastest QPCC/MPCC solver; requires installation (see nosnoc README).
model_plant_mismatch = 0; % Use a different integrator for control response simulation

%% basic timing and FESD settings
use_fesd = 1; % Use the FESD discretization method for accurate switch detection
N_FE = 2; % FESD integration steps per control interval
n_s = 2; % Number of collocation points in FESD
gamma_h = 1.0;

DT = 0.1;
T_pred           = 2.0;   % MPC prediction horizon
t_new_set_point  = 4;
t_disturbance    = 4;
T_sim            = 6.0;   % Total simulated time

% --- convert to step counts (integers) ---
N_step_new_set_point = round(t_new_set_point/DT);
N_step_disturbance   = round(t_disturbance/DT);

N_steps  = round(T_sim/DT);
N_stages = round(T_pred/DT);

tol = 1e-7;
N_homotopy = 10;


%% Define model and MPC parameters
m1 = 1.0;  % Mass of first link (kg)
m2 = 1.2;  % Mass of second link (kg)
l1 = 1.0;  % Length of first link (m)
l2 = 1.2;  % Length of second link (m)
g = 9.81;  % Gravity (m/s^2)
b1 = 0.3;  % Viscous friction coefficient at joint 1
b2 = 0.3;  % Viscous friction coefficient at joint 2
mu1 = 0.7; % Coulomb friction coefficient at joint 1
mu2 = 0.7; % Coulomb friction coefficient at joint 2

% Define states
q1 = SX.sym('q1');  % First joint angle
q2 = SX.sym('q2');  % Second joint angle
q1_dot = SX.sym('q1_dot');  % First joint angular velocity
q2_dot = SX.sym('q2_dot');  % Second joint angular velocity
x = [q1; q2; q1_dot; q2_dot];

% Define controls
tau1 = SX.sym('tau1');  % Torque at joint 1
tau2 = SX.sym('tau2');  % Torque at joint 2
u = [tau1; tau2];

% Equations of motion (derived from Euler-Lagrange)
M11 = m1*l1^2/3 + m2*(l1^2 + l2^2/3 + l1*l2*cos(q2));
M12 = m2*(l2^2/3 + l1*l2*cos(q2)/2);
M21 = M12;
M22 = m2*l2^2/3;
M = [M11, M12; M21, M22];

C1 = -m2*l1*l2*sin(q2)*(2*q1_dot*q2_dot + q2_dot^2)/2;
C2 = m2*l1*l2*sin(q2)*(q1_dot^2)/2;
C = [C1; C2];

G = [-g*(m1*l1/2 + m2*l1)*sin(q1) - g*m2*l2/2*sin(q1+q2);
    -g*m2*l2/2*sin(q1+q2)];

F_v = [-b1*q1_dot; -b2*q2_dot];  % Viscous friction

% Compute accelerations
q_ddot = M \ (u - C - G - F_v);

% Define system dynamics
f = [q1_dot; q2_dot; q_ddot];

% Define bounds for states and controls
lbx = [-2*pi; -2*pi; -inf; -inf]; % Lower bounds on states
ubx = [2*pi; 2*pi; inf; inf];     % Upper bounds on states
lbu = [-10; -10];                 % Lower bounds on controls
ubu = [10; 10];                   % Upper bounds on controls

% Define reference trajectories
x_ref = [pi; 0; 0; 0];  % Upright position
u_ref = [0; 0];         % Zero control
x0 = [0;0;0;0];

% Define weight matrices for cost function
Q = diag([10, 10, 0.1, 0.1]); % State cost weights
R = diag([0.1, 0.1]);         % Control cost weights
P = diag([100, 100, 100, 100]);


%% all options
problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
homotopy_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosnoc problems.
ccopt_options = nosnoc.ccopt.Options(); % CCopt options
mpc_options = nosnoc.mpc.Options(); % Initialize all options related to the MPC implementation.
% mpc_options.do_shift_initialization = false;
% mpc_options.warmstart_qpec = true;
mpc_options.solve_advanced_problem = 0;
mpc_options.advanced_n_qpecs = 1; % 
mpc_options.discard_constraints_in_hessian = true; % Use Gauss-Newton Hessian approximation
mpc_options.sqpec_hessian_convexification = "PROJECT";
% mpc_options.eps_hessian = 1e-9;
% mpc_options.rho_lm = 1e-3;
% Homotopy options
homotopy_options.complementarity_tol = tol; % Value used to drive the complementarity residual to zero.
homotopy_options.N_homotopy = N_homotopy; % Maximum number of homotopy iterations.
homotopy_options.print_level = 3;
% homotopy_options.sigma_0 = sigma_0;
homotopy_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; %
homotopy_options.homotopy_steering_strategy = "DIRECT";
homotopy_options.lift_complementarities = 0;

%% CCOpt settings (our fastest solver; requires installation)
ccopt_options.opts_madnlp.linear_solver = 'Ma27Solver';
ccopt_options.opts_madnlp.tol = tol;
%ccopt_options.opts_ccopt.relaxation_update.TYPE = 'ProportionalRelaxationUpdate';
%ccopt_options.opts_ccopt.relaxation_update.rolloff_slope = 2.0;
%ccopt_options.opts_ccopt.relaxation_update.rolloff_point = 1e-6;
%ccopt_options.opts_ccopt.relaxation_update.sigma_min = tol;
ccopt_options.opts_ccopt.q_regularization = 'critical_rho';
ccopt_options.opts_ccopt.critical_rho_factor = 0.9999;
%ccopt_options.opts_ccopt.q_regularization = 'eigenvalue_decomposition';
%ccopt_options.opts_ccopt.q_regularization = 'critical_rho';
%ccopt_options.opts_madnlp.barrier.TYPE = 'QualityFunctionUpdate';
%ccopt_options.opts_madnlp.barrier.mu_min = 1e-9;
%ccopt_options.opts_ccopt.use_specialized_barrier_update = true;

%% Choosing the Runge-Kutta method and number of stages
problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
problem_options.n_s = n_s; % Number of stage points in the RK method (determines accuracy)
problem_options.cross_comp_mode = "FE_FE";

% Time settings
problem_options.N_stages = N_stages; % Number of control intervals over the time horizon.
problem_options.N_finite_elements = N_FE; % Number of finite elements (integration steps) on every control interval (optionally a vector may be passed).
problem_options.T = T_pred; % Time horizon (for a time-optimal problem, the actual horizon length is obtained by solving the problem).
problem_options.use_fesd = use_fesd;
problem_options.gamma_h = gamma_h;
problem_options.euler_cost_integration = false;

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

% Equation of motion
q_ddot = M \ (u - C - G - F_v);

% Define dynamics
f_base = [q1_dot; q2_dot; q_ddot];

% f_11 = f_base + [0;0;-mu1;0];
% f_12 = f_base + [0;0;mu1;0];

f_11 = [0;0;-mu1;0];
f_12 = [0;0;mu1;0];

f_21 = [0;0;0;-mu2];
f_22 = [0;0;0;mu2];

c1 = q1_dot;
c2 = q2_dot;
% Sign matrix for the modes
S1 = [1;-1];
S2 = [1;-1];
F1 = [f_11 f_12];
F2 = [f_21 f_22];

model.f_0 = f_base; % Isolate the part of the dynamics that is the same for both modes. This reduces the nonlinearity.
model.S = {S1,S2};
model.c = {c1,c2};
model.F = {F1 F2};

model.f_q = (x-x_ref)'*Q*(x-x_ref) + (u-u_ref)'*R*(u-u_ref); % Add stage cost
model.f_q_T = (x-x_ref)'*P*(x-x_ref); % Add terminal quadratic cost

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
if model_plant_mismatch
    sim_model = nosnoc.model.Pss(); % Initialize a nosnoc model (depending on the underlying nonsmooth model, several options exist)
    % define differential states and populate the model.
    sim_model.x = x;
    sim_model.x0 = x0;
    sim_model.lbx = -inf(4,1);
    sim_model.ubx = inf(4,1);
    % define control vectors
    sim_model.u = u;
    sim_model.lbu = -inf(2,1);
    sim_model.ubu = inf(2,1);
    % Dynamics of the piecewise smooth systems
    % Define the regions of the PSS
    sim_model.f_0 = f_base;
    sim_model.c = {c1 c2};
    sim_model.S = {S1 S2};
    sim_model.F = {F1 F2}; % The columns of this matrix store the vector fields of each region.

    sim_problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
    integrator_options = nosnoc.integrator.Options();
    integrator_options.print_level = 1;
    sim_solver_options = integrator_options.fesd_solver_opts; % Initialize all options related to the MPEC solver used for solving nosnoc problems.
    % Choose the Runge-Kutta method and number of stages
    sim_problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
    sim_problem_options.n_s = n_s; % Number of stage points in the RK method (determines accuracy)

    % Time settings
    sim_problem_options.N_finite_elements = 2; % Number of finite elements (integration steps) on every control interval (optionally a vector may be passed).
    sim_problem_options.T_sim = problem_options.h;
    sim_problem_options.N_sim = 5;
    sim_problem_options.print_level = 0;

    % Simulation solver options
    sim_solver_options.print_level = 0;
    sim_solver_options.homotopy_steering_strategy = HomotopySteeringStrategy.ELL_INF; % Use the $\ell_{\infty}$ steering strategy
    sim_solver_options.complementarity_tol = 1e-10; % Value used to drive the complementarity residual to zero.
    % sim_solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; % Using an HSL solver is a significant speedup over the default 'mumps', but requires installation
    integrator = nosnoc.Integrator(sim_model, sim_problem_options, integrator_options);
end
%% MPC Loop

plot_intermediate_solutions = true;
if plot_intermediate_solutions
    q_plot = subplot(311); hold on;
    xlim([0 (N_steps+N_stages)*problem_options.h])
    ylim([-pi*1.1 pi*1.1])
    xlabel(q_plot, '$t$')
    ylabel(q_plot, '$q$')

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
qp_solved = []; qp_accepted = []; qpec_solved = [];
feedback_times = []; preparation_times = [];
f_opt = []; % Value function


for step=1:N_steps
    stats_prep = mpc.do_preparation();
    tp_i = stats_prep.preparation_time;
    % What would the plant do? Compute plant-model mismatch for the current MPC input
    preparation_times = [preparation_times, tp_i];

    [u_i, stats_fb] = mpc.get_feedback(x0);
    % f_opt = [f_opt, mpc.get_objective()];
    tf_i = stats_fb.feedback_time;
    feedback_times = [feedback_times, tf_i];
    if use_rtmpc
        if save_qpecs && mod(step-1, 5) == 0
            save_qpec(mpc.qpec, ['ACROBOT_CONVEX_', num2str(step)])
        end
        qp_solved = [qp_solved, stats_fb.qp_solved];
        qp_accepted = [qp_accepted, stats_fb.qp_accepted];
        qpec_solved = [qpec_solved, stats_fb.qpec_solved];
    else
        qp_solved = [qp_solved, nan];
        qp_accepted = [qp_accepted, nan];
        qpec_solved = [qpec_solved, nan];
    end


    fprintf("MPC step: %d/%d, Preparation time: %d [s], Feedback time: %d [s]\n", step, N_steps, tp_i, tf_i);
    fprintf("Control input: %2.2f\n", u_i);

    f_opt = [f_opt, mpc.get_objective()];

    if model_plant_mismatch
        integrator.set_x0(x0);
        [t_grid, x_sim] = integrator.simulate("u", repmat(u_i, [1,5]), "x0", x0);
        x0 = x_sim(:, end);
    else
        x0 = mpc.get_predicted_state();
        % if mod(N_steps,6) == 0
        %  x0(1:2) = x0(1:2) + 0.1-0.2*rand(2,1);
        % end
    end


    if step == N_step_disturbance
        % Add a large disturbance to the angle
        x0(2) = x0(2)-1;
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

        cla(q_plot); hold on;
        cla(v_plot); hold on;
        cla(u_plot); hold on;

        plot(q_plot, t, x(1:2,:))
        plot(q_plot, t(end)+t_grid, q_res)
        yline(q_plot, x_ref(1:2), 'k--')
        xlabel(q_plot, '$t$')
        ylabel(q_plot, '$q$')

        plot(v_plot, t, x(3:4,:))
        plot(v_plot, t(end)+t_grid, v_res)
        yline(v_plot, x_ref(3:4), 'k--')
        xlabel(v_plot, '$t$')
        ylabel(v_plot, '$v$')

        if ~isempty(u)
            u_hist = [u, u(:,end)];
            stairs(u_plot, t, u_hist')
        end
        stairs(u_plot, t(end)+t_grid_u, [u_res, u_res(:,end)]')
        yline(u_plot, u_ref, 'k--')
        xlabel(u_plot, '$t$')
        ylabel(u_plot, '$u$')
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


%% Visualization

% end
x_sim = x;
N = length(x);
% Visualization with video recording
save_video  = 1;
if save_video
    v = VideoWriter('acrobot.mp4', 'MPEG-4');
    % v = VideoWriter('acrobot.avi');
    % v.FrameRate = 1;
    open(v);
end
figure; hold on;
axis([-2.5 2.5 -2.5 2.5]);
grid on;
for k = 1:N
    clf; hold on;
    theta1 = x_sim(1,k);
    theta2 = x_sim(2,k);
    p1 = [l1*sin(theta1); -l1*cos(theta1)];
    p2 = p1 + [l2*sin(theta1 + theta2); -l2*cos(theta1 + theta2)];
    plot([0 p1(1)], [0 p1(2)], 'r', 'LineWidth', 2);
    plot([p1(1) p2(1)], [p1(2) p2(2)], 'b', 'LineWidth', 2);
    plot(p1(1), p1(2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    plot(p2(1), p2(2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    xlim([-2.5 2.5]); 
    ylim([-2.5 2.5]);
    axis equal;
    if save_video
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
    pause(0.01);
end
if save_video
    close(v);
end
%% check:
N_check = N_steps-1;
fprintf(['QP solved: %.1f%%  |  QP accepted: %.1f%%  |  QPEC solved (fallback): %.1f%%\n'], 100*sum(qp_solved)/N_check, 100*sum(qp_accepted)/N_check, 100*sum(qpec_solved)/N_check);

%%
figure
stairs(preparation_times,'DisplayName','Preparation time','LineWidth',2)
hold on
stairs(feedback_times,'DisplayName','Feedback time','LineWidth',2)
legend show