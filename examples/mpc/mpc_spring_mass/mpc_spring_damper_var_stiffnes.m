%% Minimal MPC example for a switched spring-mass system.
% This script sets up and solves an MPC problem for a piecewise-smooth
% oscillator and compares real-time MPC with a standard full MPC variant.
% Here, mpecopt is used as qpcc solver.

clear; clc; close all;
import casadi.*
import nosnoc.*
%%
model_plant_mismatch = 0;
use_rtmpc = 1;
%% Parameters
N_stages = 20; % Number of control intervals
N_steps = 20; % Number of closed-loop MPC steps


%% Model and MPC parameters
m = 1.0;       % Mass (kg)
k1 = 0.5;      % Spring constant in the negative region (N/m)
k2 = 3.0;      % Spring constant in the positive region (N/m)
c = 0.7;       % Damping coefficient (Ns/m)

c_v = 0.005; % Coefficient of air resistance
c_v  = 0;

Q = diag([1.0, 0.5]);  % State tracking weight
R = 0.1;               % Control input weight
P = diag([5.0, 2.5]);  % Terminal cost weight

x1_max = 3.0;
x2_max = 4.0;
% Input constraints
u_max = 1.0;
x0 = [-1.0; 0.0];
x_ref = [1.5; 0.0];    % Reference state [position; velocity]
% x_ref = [0; 0.0];    % Reference state [position; velocity]

u_ss1 = -k1/m*x_ref(1) - c/m*x_ref(2);
u_ss2 = -k2/m*x_ref(1) - c/m*x_ref(2);
u_max = 5; % Set point reachable
% u_max = 2; % Set point unreachable
u_ref = -1*u_ss1;
u_ref = -u_ss2;
% u_ref = 0;

%% all options
problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
homotopy_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosnoc problems.
mpecopt_options = mpecopt.Options(); % Initialize all options related to the MPEC solver used for solving nosnoc problems.
mpc_options = nosnoc.mpc.Options(); % Initialize all options related to the MPC implementation.

%%
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
model.lbx = [-x1_max;-x2_max]; % Lower bounds on states
model.ubx = [x1_max;x2_max]; % Upper bounds on states
% define control vectors
u = SX.sym('u');  % CasADi symbolic variables for controls
model.u = u;
model.lbu = -u_max;
model.ubu = u_max;

rhs_neg = [x(2); 
           -k1/m*x(1) - c/m*x(2) + 1/m*u];
rhs_pos = [x(2); 
           -k2/m*x(1) - c/m*x(2) + 1/m*u];

model.c = x(1); 
model.S = [-1;1]; % Region R_1 is defined by c<0 (hence the -1), region R_2 by c>0 (hence the +1) in the sign matrix S.
model.F = [rhs_neg rhs_pos]; % The columns of this matrix store the vector fields of each region.


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
% solver_options.homotopy_steering_strategy ="ELL_INF";
homotopy_options.lift_complementarities = 0;

%% Optimization solver
mpecopt_options.compute_tnlp_stationary_point = false;
mpecopt_options.settings_lpec.lpec_solver = "Gurobi";
mpecopt_options.settings_casadi_nlp.ipopt.linear_solver = 'ma27';
mpecopt_options.relax_and_project_homotopy_parameter_steering = "Direct";
mpecopt_options.rho_TR_phase_i_init = 1e-2;
mpecopt_options.rho_TR_phase_ii_init = 1e-5;
mpecopt_options.verbose_summary = 0;

%% some mpc options
mpc_options.fast_sigma_0 = 1e-1;
mpc_options.do_shift_initialization = true;
mpc_options.solve_advanced_problem = false;

%% create mpc object
if use_rtmpc
    mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, mpecopt_options, mpecopt_options);
else
    mpc = nosnoc.mpc.FullMpcc(model, mpc_options, problem_options, homotopy_options);
end

%% Create sim model and integrator
% Create the true model with air resistance
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


% Mode 1: Negative position (x1 < 0) with k1 stiffness
rhs_neg = [x(2); 
           -k1/m*x(1) - c/m*x(2) + 1/m*u+ c_v*x(2)^2];
                
% Mode 2: Positive position (x1 >= 0) with k2 stiffness
rhs_pos = [x(2); 
           -k2/m*x(1) - c/m*x(2) + 1/m*u+ c_v*x(2)^2];

% Define the regions of the PSS
sim_model.c = x(1);
sim_model.S = [-1;1]; % Region R_1 is defined by c<0 (hence the -1), region R_2 by c>0 (hence the +1) in the sign matrix S.
sim_model.F = [rhs_neg rhs_pos]; % The columns of this matrix store the vector fields of each region.

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
%% Do MPC with accurate model

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