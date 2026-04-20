%% Minimal MPC example for a spring-mass system with Coulomb friction.
% The MPC problem is to drive the mass from its initial state to a desired
% position and velocity while respecting input and state bounds and keeping
% the control effort small. The model is piecewise smooth because the
% friction force changes sign with the velocity, which creates a switch at
% zero velocity.

clear; clc; close all;
import casadi.*
import nosnoc.*
%%
use_rtmpc = 1;
model_plant_mismatch = 1;
%% Parameters
N_stages = 20;
N_steps = 40; % Number of closed-loop MPC steps

m = 1.0;      % Mass (kg)
k = 0.8;      % Spring constant (N/m)
c = 0.5;      % Damping coefficient (Ns/m)
Fc = 0.2;     % Coulomb friction coefficient (N)

c_v = 0.005; % Coefficient of air resistance
c_v  = 0;

Q = diag([1.0, 1.0]);  % State tracking weight
R = 1e-1;              % Control input weight
P = Q;                 % Terminal cost weight

x1_min = -3.0;
x1_max = 3.0;
% Input constraints
u_max = 1.0;

x0 = [0;0]; % Initial value
x_ref = [2; 0.0];    % Reference state [position; velocity]
% x_ref = [1.25;0];

u_ss = -k/m*x_ref(1) - c/m*x_ref(2) - Fc/m;
u_max = 2; % Set point reachable
% u_max = 1; % Set point unreachable
u_ref = -1*u_ss;
% u_ref = 0;

%% all options
problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
homotopy_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosnoc problems.
mpecopt_options = mpecopt.Options(); % Initialize all options related to the MPEC solver used for solving nosnoc problems.
mpc_options = nosnoc.mpc.Options(); % Initialize all options related to the MPC implementation.

%% problem definition
% Choose the Runge-Kutta method and number of stages
problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
problem_options.n_s = 2; % Number of stage points in the RK method (determines accuracy)
problem_options.cross_comp_mode = "FE_FE";

% Time settings
problem_options.N_stages = N_stages; % Number of control intervals over the time horizon.
problem_options.N_finite_elements = 2; % Number of finite elements (integration steps) on every control interval (optionally a vector may be passed).
problem_options.T = 3;    % Time horizon (for a time-optimal problem, the actual horizon length is obtained by solving the problem).

model = nosnoc.model.Pss(); % Initialize a nosnoc model (depending on the underlying nonsmooth model, several options exist)


%% model

% define differential states and populate the model.
x = SX.sym('x',2);
model.x = x;
model.x0 = x0;
v_max = 5;
model.lbx = [x1_min;-v_max]; % Lower bounds on states
model.ubx = [x1_max;v_max]; % Upper bounds on states
% define control vectors
u = SX.sym('u');  % CasADi symbolic variables for controls
model.u = u;
model.lbu = -u_max;
model.ubu = u_max;

% Mode 1: Positive velocity (x2 > 0)
rhs_positive = [x(2);
    -k/m*x(1) - c/m*x(2) - Fc/m*1 + 1/m*u];

% Mode 2: Negative velocity (x2 < 0)
rhs_negative = [x(2);
    -k/m*x(1) - c/m*x(2) - Fc/m*(-1) + 1/m*u];

model.c = x(2); % Single switching function (implies two regions)
model.S = [-1;1]; % Region R_1 is defined by c<0 (hence the -1), region R_2 by c>0 (hence the +1) in the sign matrix S.
model.F = [rhs_negative rhs_positive]; % The columns of this matrix store the vector fields of each region.


model.f_q = (x-x_ref)'*Q*(x-x_ref) + R*(u-u_ref)^2; % Add stage cost
model.f_q_T = (x-x_ref)'*P*(x-x_ref); % Add terminal quadratic cost

% Setup solver an mpc options
% solver_options.homotopy_update_rule = 'superlinear'; % Use a superlinear update rule for the relaxation parameter sigma.
% solver_options.homotopy_update_slope = 0.05; % Rate at which the relaxation sigma is reduced: sigma_i+1 = kappa*sigma_i
% solver_options.homotopy_update_exponent = 2; % Rate at which the relaxation sigma is reduced: sigma_i+1 = sigma_i^kappa
homotopy_options.complementarity_tol = 1e-6;
homotopy_options.N_homotopy = 10;
%solver_options.print_level = 0;
homotopy_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; % Using an HSL solver is a significant speedup over the default 'mumps', but requires installation
% homotopy_options.opts_casadi_nlp.ipopt.hsllib = '/home/anton/tools/HSL_jll.jl-2023.11.7/override/lib/x86_64-linux-gnu-libgfortran5/libhsl.so';
% homotopy_options.homotopy_steering_strategy ="ELL_INF";
homotopy_options.lift_complementarities = 0;
%%  mpc opts
mpc_options.solve_advanced_problem = false;
mpc_options.fast_sigma_0 = 1e-5;
mpc_options.do_shift_initialization = true;
%% QPCC/QPEC solver options
qpcc_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosnoc problems.
% qpcc_options.sigma_0 = 1e-2;
% qpcc_options.homotopy_update_rule = 'superlinear';
qpcc_options.homotopy_update_slope = 0.01;
qpcc_options.homotopy_steering_strategy = "DIRECT";
% qpcc_options.opts_casadi_nlp.ipopt.mu_strategy = 'monotone';
qpcc_options.opts_casadi_nlp.ipopt.tol = 1e-9;

% lcqpow_options = LCQPow_options(); % TODO(@anton) please make a class for this.
% lcqpow_options.printLevel = 1;
% lcqpow_options.qpSolver = 1;
% lcqpow_options.maxIterations = 1000;
% lcqpow_options.initialPenaltyParameter = 10;
% lcqpow_options.penaltyUpdateFactor = 10;
% lcqpow_options.complementarityTolerance = 1e-6;
% lcqpow_options.stationarityTolerance = 1e-6;
%lcqpow_options.qpOASES_options.printLevel = ;

gurobi_options = nosnoc.qpec.GurobiOptions();
% gurobi_options.method = 'reg';
% gurobi_options.method = 'miqp';
gurobi_options.method = 'sos1';


%% create mpc object
if use_rtmpc
    mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, homotopy_options,homotopy_options);
else
    mpc = nosnoc.mpc.FullMpcc(model, mpc_options, problem_options, homotopy_options);
end

%% Create sim model and integrator
sim_model = nosnoc.model.Pss(); % Initialize a nosnoc model (depending on the underlying nonsmooth model, several options exist)

% define differential states and populate the model.
x = SX.sym('x',2);
sim_model.x = x;
sim_model.x0 = x0;
sim_model.lbx = -inf(2,1);
sim_model.ubx = inf(2,1);
% define control vectors
u = SX.sym('u');  % CasADi symbolic variables for controls
sim_model.u = u;
sim_model.lbu = -inf;
sim_model.ubu = inf;
% Dynamics of the piecewise smooth systems

rhs_positive = [x(2);
    -k/m*x(1) - c/m*x(2) - Fc/m*1 + 1/m*u - c_v*x(2)^2];
% Mode 2: Negative velocity (x2 < 0)
rhs_negative = [x(2);
    -k/m*x(1) - c/m*x(2) - Fc/m*(-1) + 1/m*u + c_v*x(2)^2];

% Define the regions of the PSS
sim_model.c = x(2);
sim_model.S = [-1;1]; % Region R_1 is defined by c<0 (hence the -1), region R_2 by c>0 (hence the +1) in the sign matrix S.
sim_model.F = [rhs_negative rhs_positive]; % The columns of this matrix store the vector fields of each region.

sim_problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
integrator_options = nosnoc.integrator.Options();
sim_solver_options = integrator_options.fesd_solver_opts; % Initialize all options related to the MPEC solver used for solving nosnoc problems.

% Choose the Runge-Kutta method and number of stages
sim_problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
sim_problem_options.n_s = 4; % Number of stage points in the RK method (determines accuracy)

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

%% MPC Loop
plot_intermediate_solutions = true;
if plot_intermediate_solutions
    q_plot = subplot(311); hold on;
    xlim([0 (N_steps+N_stages)*problem_options.h])
    ylim([-1 2])
    xlabel(q_plot, '$t$')
    ylabel(q_plot, '$q$')

    v_plot = subplot(312); hold on;
    xlim([0 (N_steps+N_stages)*problem_options.h])
    ylim([-1.2 1.2])
    xlabel(v_plot, '$t$')
    ylabel(v_plot, '$v$')

    u_plot = subplot(313); hold on;
    xlim([0 (N_steps+N_stages)*problem_options.h])
    ylim([-1.05*u_max 1.05*u_max])
    xlabel(u_plot, '$t$')
    ylabel(u_plot, '$u$')
end

x = model.x0; u = []; t = 0; tf = [];
x0 = x;
for step=1:N_steps
    [u_i, stats] = mpc.get_feedback(x0);
    tf_i = stats.feedback_time;
    tf = [tf, tf_i];
    fprintf("MPC step: %d, Feedback time: %d\n", step, tf_i);
    mpc.do_preparation();
    if model_plant_mismatch
        integrator.set_x0(x0);
        [t_grid, x_sim] = integrator.simulate("u", repmat(u_i, [1,5]), "x0", x0);
        x0 = x_sim(:, end);
    else
        x0 = mpc.get_predicted_state();
    end

    % Plot intermediate solution by getting x, t, and u from mpc object
    % Using `.get`
    if plot_intermediate_solutions
        x_res = mpc.get('x');
        q_res = x_res(1,:);
        v_res = x_res(2,:);
        t_grid = mpc.get_time_grid();
        u_res = mpc.get('u');
        t_grid_u = mpc.get_control_grid();
        %
        cla(q_plot); hold on;
        cla(v_plot); hold on;
        cla(u_plot); hold on;

        plot(q_plot, t, x(1,:))
        plot(q_plot, t(end)+t_grid, q_res)
        yline(q_plot,x_ref(1),'k--')
        xlabel(q_plot, '$t$')
        ylabel(q_plot, '$q$')

        plot(v_plot, t, x(2,:))
        plot(v_plot, t(end)+t_grid, v_res)
        yline(v_plot,x_ref(2),'k--')
        xlabel(v_plot, '$t$')
        ylabel(v_plot, '$v$')

        if step>1
            stairs(u_plot, t, [u, u(end)])
        end
        stairs(u_plot, t_grid_u+t(end), [u_res, u_res(end)])
        yline(u_plot,u_ref,'k--')
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
plot(t,x(1,:),'LineWidth',1.5)
hold on;
yline(x_ref(1),'k--','LineWidth',1.5)
xlabel("$t$")
ylabel("$q$")
subplot(312)
plot(t,x(2,:),'LineWidth',1.5)
hold on;
yline(x_ref(2),'k--','LineWidth',1.5)
xlabel("$t$")
ylabel("$v$")
% ylim([-5 25])
subplot(313)
stairs(t,[u,u(end)],'LineWidth',1.5)
hold on
yline(u_max,'k--','LineWidth',1.5)
yline(-u_max,'k--','LineWidth',1.5)
xlabel("$t$")
ylabel("$u$")
ylim([-1.1*u_max 1.1*u_max])

% figure
% bar(tf)
% ylabel("feedback time (s)")
% xlabel("step")