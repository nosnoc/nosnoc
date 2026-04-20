clear; clc; close all;
import casadi.*
%% Description
% MPC for an ideal gas-liquid tank system with a pressure control valve.
% The MPC problem regulates the valve opening to control the tank behavior
% while penalizing valve usage and liquid holdup. The model is piecewise
% smooth because the outlet dynamics switch depending on whether the liquid
% level is above or below the separator threshold.
%
% Example from: Globally Convergent Method for Optimal Control of Hybrid
% Dynamical Systems by Saif R. Kazi, Mandar Thombre, Lorenz Biegler
% 12th IFAC International Symposium on Advanced Control of Chemical
% Processes, July 14-17, 2024, Toronto, Canada.

%% plot and video settings
latexify_plot()
nice_plot_colors;

%%
filename_pdf = 'liquid_gas_mpc';
linewidth = 1.5;
use_rtmpc = 1; % Use real-time MPC instead of solving each OCP to full convergence
use_ccopt = 0; % Recommended: use our fastest QPCC/MPCC solver; requires installation (see nosnoc README)
save_qpecs = 0; % Save QPEC/QPCC subproblems

%% Model settings
% Discretization and simulation settings
model_plant_mismatch = 1; % Use a separate plant integrator for closed-loop simulation
model_plant_same_integrator = 0; % If true, use the same discretization/integrator settings for plant and MPC model

N_sim = 10;
N_steps = 50; % Number of MPC steps in simulation


N_stages = 50; % Number of control stages in the prediction horizon
DT = 0.25;

N_FE = 1; % Number of finite elements per control interval
use_fesd = 1; % Use FESD discretization
gamma_h = 0.0; % Step-size regularization parameter for FESD

%% solver options
tol = 1e-6; % Complementarity tolerance
N_homotopy = 6; % Maximum number of homotopy iterations
sigma_0 = 1; % Initial homotopy relaxation parameter


%% nosnoc setup
problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
mpc_options = nosnoc.mpc.Options(); % Initialize all options related to the MPC implementation.
homotopy_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosnoc problems.

ccopt_options = nosnoc.ccopt.Options(); % CCopt options
%ccopt_options.opts_ccopt.relaxation_update.TYPE = 'RelaxLBUpdate';
%ccopt_options.opts_ccopt.relaxation_update.TYPE = 'RolloffRelaxationUpdate';
ccopt_options.opts_madnlp.linear_solver = 'Ma27Solver';
ccopt_options.opts_madnlp.tol = tol;

%% Model variables and parameters
% Specify initial and desired state

% CasADi symbolic variables
n_x = 2;
n_u = 1;
M_G = SX.sym('M_G'); % Vapor holdup
M_L = SX.sym('M_L'); % Liquid holdup

x = [M_G;M_L];
u = SX.sym('u',n_u); % Control
% Parameters
F_L = 2.5; % mol/sec
F_G = 0.1; % mol/sec
V = 10;
Vs = 5; % liters
T = 300; % K
P_out = 1; % atm
rho_L = 50;
k_L = 1;
k_G = 1;
% Initial values
L0 = 0.25;
G0 = 0.25;
P0 = 35;
M_L0 = 260;
M_G0 = 6.83;
% R = 8.314; % J/mol K
R = 0.0821; % L atm/(mol K)

%% Model

c = M_L/rho_L-Vs;

% Ideal gas equation
% M_G*R*T/P + M_L/rho_L = V;
% M_G*R*T/P = V - M_L/rho_L;
% P/(M_G*R*T) = 1/(V - M_L/rho_L);
P = M_G*R*T/(V-M_L/rho_L);

P_fun = Function('P_fun', {x}, {P});
c_fun = Function('P_fun', {x}, {c});

L = k_L*u*(P-P_out);
G = k_G*u*(P-P_out);

% c > 0: liquid outlet active
f_1 = [F_G;F_L-L];

% c < 0: gas outlet active
f_2 = [F_G-G; F_L];

ubx = [inf; inf];
lbx = [-inf; -inf];
u_max = 0.15;
u_ref = 0.1;

% Stage cost
f_q = 100*(u-u_ref)^2 + M_L;
% Terminal cost
f_q_T  = 10*M_L;
F = [f_1, f_2];
x0 = [M_G0;M_L0];
% Populate nosnoc model
model = nosnoc.model.Pss();
model.c = c;         % Switching function
model.S = [1; -1];   % Sign matrix: f_1 for c>0, f_2 for c<0
model.F = F;
model.f_q = f_q;
model.f_q_T = f_q_T;
model.lbx = lbx;
model.ubx = ubx;
model.x = x;
model.x0 = x0;
model.u = u;
model.lbu = 0;
model.ubu = u_max;

% model.g_terminal = model.x-x_ref;
% problem_options.relax_terminal_constraint = "ELL_INF";

%% OCP problem options
problem_options.T = N_stages*DT; % Time horizon
problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Runge-Kutta scheme
problem_options.n_s = 2; % Number of RK stage points
problem_options.dcs_mode = 'Stewart'; % DCS formulation
problem_options.N_stages = N_stages; % Number of control intervals
problem_options.N_finite_elements = N_FE; % Number of finite elements per control interval
problem_options.cross_comp_mode = "FE_FE"; % Cross complementarity handling mode
problem_options.use_fesd = use_fesd; % Enable FESD discretization
problem_options.gamma_h = gamma_h; % Step-size regularization parameter

%% MPEC solver options
homotopy_options.complementarity_tol = tol; % Value to drive the complementarity residual to.
homotopy_options.N_homotopy = N_homotopy; % Maximum number of homotopy iterations.
homotopy_options.print_level = 3; % Solver verbosity
homotopy_options.sigma_0 = sigma_0; % Initial homotopy relaxation parameter
homotopy_options.opts_casadi_nlp.ipopt.linear_solver = 'mumps'; %
homotopy_options.homotopy_steering_strategy = "ELL_INF"; % Homotopy steering strategy

%% MPC options
mpc_options.fast_sigma_0 = 1e-1; % Initial sigma used for fast MPC warm starts
mpc_options.do_shift_initialization = true; % Shift previous solution for initialization
mpc_options.solve_advanced_problem = false; % Disable advanced auxiliary MPC problem
mpc_options.sqpec_hessian_convexification = "MIRROR"; % Hessian convexification strategy

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
%% create plant
if model_plant_mismatch
    % Integrator model
    sim_model = nosnoc.model.Pss();
    sim_model.c = c;
    sim_model.x = x;
    sim_model.x0 =  x0;
    sim_model.u = u;
    sim_model.F = F;
    sim_model.S = [1;-1];
    sim_model.lbx = -inf(n_x,1); sim_model.ubx = inf(n_x,1);
    sim_model.lbu = -inf;  sim_model.ubu = inf;
    sim_problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
    integrator_options = nosnoc.integrator.Options();
    integrator_options.integrator_plugin = "FESD"; % Use the FESD integrator

    sim_solver_options = integrator_options.fesd_solver_opts; % The FESD integrator uses an MPEC solver; modify its options here
    % Choosing the Runge-Kutta method and number of stages
    sim_problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
    if model_plant_same_integrator
        sim_problem_options.n_s = 2; % Number of stage points in the RK method (determines accuracy)
        sim_problem_options.N_finite_elements = N_FE; % Number of finite elements (integration steps) on every control interval (optionally a vector might be passed).
        sim_problem_options.N_sim = 1; % Number of simulation intervals per sample
        N_sim = 1;
    else
        sim_problem_options.n_s = 3; % Number of stage points in the RK method (determines accuracy)
        sim_problem_options.N_finite_elements = 2; % Number of finite elements (integration steps) on every control interval (optionally a vector might be passed).
        sim_problem_options.N_sim = N_sim; % Number of simulation intervals per sample
    end
    sim_problem_options.T_sim = DT; % Simulation horizon per MPC step
    sim_solver_options.print_level = 0; % Solver verbosity
    sim_solver_options.homotopy_steering_strategy = HomotopySteeringStrategy.ELL_INF; % Use the $\ell_{\infty}$ steering strategy
    sim_solver_options.complementarity_tol = 1e-7; % Value to drive the complementarity residual to.
    sim_solver_options.opts_casadi_nlp.ipopt.linear_solver = 'mumps'; % Using an HSL solver is a significant speedup over the default 'mumps' but requires installation

    % create nosnoc integrator
    integrator = nosnoc.Integrator(sim_model, sim_problem_options, integrator_options);
end

%% MPC loop
x0 = model.x0;
u = [];
t = 0;
tf = [];
x = x0;
d_hat = zeros(4,1);
% error_in_x0 = 0;
error_in_x0 = zeros(4,1);
error_in_predicted_states = [];
f_opt = [];
x_high_res = x0;
t_high_res = 0;
x_mpc_pred = x0;

plot_intermediate_solutions = 1;
if plot_intermediate_solutions
    M_G_plot = subplot(311); hold on;
    xlim([0 (N_steps+N_stages)*problem_options.h])
    xlabel(M_G_plot, '$t$')
    ylabel(M_G_plot, '$M_G$')

    M_L_plot = subplot(312); hold on;
    xlim([0 (N_steps+N_stages)*problem_options.h])
    xlabel(M_L_plot, '$t$')
    ylabel(M_L_plot, '$P$')

    u_plot = subplot(313); hold on;
    xlim([0 (N_steps+N_stages)*problem_options.h])
    ylim([0.3*u_max u_max])
    xlabel(u_plot, '$t$')
    ylabel(u_plot, '$u$')
end

for step=1:N_steps
    [u_i, stats] = mpc.get_feedback(x0);
    % f_opt = [f_opt, mpc.get_objective()];
    tf_i = stats.feedback_time;
    tf = [tf, tf_i];
    if save_qpecs && mod(step-1, 5) == 0
        save_qpec(mpc.qpec, ['GAS_LIQUID_CONVEX_', num2str(step)])
    end
    fprintf("MPC step: %d, Feedback time: %d\n", step, tf_i);
    fprintf("Control input: %2.2f\n", u_i);
    mpc.do_preparation();
    % What would the plant do? Compute plant-model mismatch for the current MPC input

    if model_plant_mismatch
        % Advance plant in time
        integrator.set_x0(x0);
        [t_grid_sim, x_sim] = integrator.simulate("u", repmat(u_i, [1,N_sim]), "x0", x0);
        x_high_res = [x_high_res,x_sim];
        t_high_res = [t_high_res, t_high_res(end) + t_grid_sim];
        x0 = x_sim(:, end);
    else
        x0 = mpc.get_predicted_state();
    end

    if plot_intermediate_solutions
        x_res = mpc.get('x');
        M_G_res = x_res(1,:);
        M_L_res = x_res(2,:);
        t_grid = mpc.get_time_grid();
        u_res = mpc.get('u');
        t_grid_u = mpc.get_control_grid();
        %
        cla(M_G_plot); hold on;
        cla(M_L_plot); hold on;
        cla(u_plot); hold on;

        plot(M_G_plot, t, x(1,:))
        plot(M_G_plot, t(end)+t_grid, M_G_res)
        xlabel(M_G_plot, '$t$')
        ylabel(M_G_plot, '$M_G$')
        
        P = full(P_fun(x));
        P_res = full(P_fun(x_res));
        % plot(M_L_plot, t, x(2,:))
        % plot(M_L_plot, t(end)+t_grid, M_L_res)
        plot(M_L_plot, t, P)
        plot(M_L_plot, t(end)+t_grid, P_res)
        xlabel(M_L_plot, '$t$')
        ylabel(M_L_plot, '$P$')
        % yline(M_L_plot,x_ref(3:4),'k--')

        if step>1
            stairs(u_plot, t, [u, u(:,end)]')
        end
        stairs(u_plot, t_grid_u+t(end), [u_res, u_res(:,end)]')
        yline(u_plot,u_ref,'k--')
        xlabel(u_plot, '$t$')
        ylabel(u_plot, '$u$')

    end

    x = [x, x0];
    u = [u,u_i];
    t = [t, t(end) + problem_options.h];
end

%% Plot
P = full(P_fun(x));
c = full(c_fun(x));
f = figure;
subplot(411)
plot(t,x(1,:),'LineWidth',linewidth)
xlabel("$t$")
ylabel("$M_G$")
grid on
yline(0,'k--','LineWidth',1.5)
subplot(412)
plot(t,x(2,:),'LineWidth',linewidth )
xlabel("$t$")
ylabel("$M_L$")
grid on
subplot(413)
plot(t,P,'LineWidth',linewidth )
xlabel("$t$")
ylabel("$P$")
grid on

subplot(414)
stairs(t,[u,u(end)],'LineWidth',linewidth)
xlabel("$t$")
ylabel("$u$")
grid on
ylim([0.07 0.12])
hold on

figure
plot(t,c,'LineWidth',linewidth)
xlabel("$t$")
ylabel("$c(x(t))$")
grid on