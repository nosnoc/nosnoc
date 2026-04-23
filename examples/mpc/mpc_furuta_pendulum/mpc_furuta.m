%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adapted for nosnoc / original acados model by Jörg Fischer
% https://github.com/Jo-Fischer/acados-STM32-NUCLEO-H745ZI/tree/master
% Description:
%   Closed-loop MPC for the Furuta Pendulum swing-up using nosnoc.
%   Coulomb friction is treated correctly as a PSS 
%
%   States:  x = [theta1, theta1_dot, theta2, theta2_dot]
%     theta1     : actuator arm angle [rad]   (horizontal plane)
%     theta1_dot : arm angular velocity [rad/s]
%     theta2     : pendulum angle [rad]  (0 = hanging DOWN, pi = upright)
%     theta2_dot : pendulum angular velocity [rad/s]
%
%   Input:   u  motor voltage [V] in [-U_max, U_max]
%
%   Cost — exact replication of acados NONLINEAR_LS cost:
%     output  y  = [pi*(1+cos(theta1/2)); theta1p; pi*(1+cos(theta2/2)); theta2p; u]
%     yref       = [pi; 0; pi; 0; 0]
%     => err_vec = y - yref = [pi*cos(theta1/2); theta1p; pi*cos(theta2/2); theta2p]
%     stage cost = err_vec' * W_x * err_vec + W_u * u^2
%     Verify: theta=0 (hanging) -> pi*cos(0)   = pi  -> max cost
%             theta=pi (upright) -> pi*cos(pi/2) = 0   -> zero cost ✓
%
%   PSS — additive Coulomb correction 
%     model.f_0 = f_base  (smooth dynamics, no Coulomb)
%     F{1}, F{2} = purely the signed Coulomb acceleration increments
%     Switch 1: sign(theta1p) — arm joint Coulomb (b1coul)
%     Switch 2: sign(theta2p) — pendulum joint Coulomb (b2coul)
%
%   Note on Furuta vs acrobot Coulomb structure:
%     In the acrobot, mu_i enters directly as a generalised force on ddq_i.
%     In the Furuta pendulum, b_coul enters the RHS *before* M^{-1}, so
%     the state-space increment is M\[±b_coul; 0] or M\[0; ±b_coul].
%     Because M is off-diagonal (coupling), each Coulomb force affects
%     BOTH accelerations. The increment is symbolic in theta2.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;

import casadi.*
import nosnoc.*
latexify_plot();

%% -----------------------------------------------------------------------
%  Flags
% ------------------------------------------------------------------------
use_rtmpc            = 1;    % 1 = Rtmpc (SQPCC-HyRTI), 0 = FullMpcc
use_ccopt            = 0;    % Recommended for best perfomance but needs CCOpt installation
model_plant_mismatch = 0;    % 1 = use separate nosnoc integrator as plant
save_video           = 1;    % save animation as furuta_mpc.mp4
plot_intermediate    = 1;    % live plot during MPC loop

%% -----------------------------------------------------------------------
%  Physical Parameters  (from Parameters_FurutaPendulum_Coulomb.m)
% ------------------------------------------------------------------------
grav     = 9.81;
kSteller = 1.0983;   % PWM unit voltage factor
kM       = 0.0236;   % motor torque constant [Nm/A]
RM       = 0.5200;   % motor resistance [Ohm]
b1vis    = 6.0700e-04;  % arm viscous friction [Nm*s/rad]
b1coul   = 0.0104;      % arm Coulomb friction [Nm]
hJ0      = 0.0021;      % arm + hanging pendulum inertia about motor axis [kg*m^2]
L1       = 0.1215;      % arm length [m]
m2       = 0.0461;      % pendulum mass [kg]
l2       = 0.0790;      % pendulum CoM distance from joint [m]
b2vis    = 2.6830e-05;  % pendulum viscous friction [Nm*s/rad]
b2coul   = 7.5418e-04;  % pendulum Coulomb friction [Nm]
hJ2      = 3.7054e-04;  % pendulum inertia about joint [kg*m^2]

%% -----------------------------------------------------------------------
%  MPC Timing & Weights
% ------------------------------------------------------------------------
DT       = 0.04;             % control interval [s]
T_pred   = 0.4;              % prediction horizon [s]
T_sim    = 2;              % total simulation time [s]

N_stages = round(T_pred/DT); % horizon steps
N_steps  = round(T_sim /DT); % MPC loop iterations

N_FE     = 2;     % FESD finite elements per interval
n_s      = 2;     % Radau IIA stages
use_fesd = 1;
gamma_h  = 1.0;

tol        = 1e-6;
N_homotopy = 10;

U_max  = 10.0;               % [V] input constraint
x0_val = [0.1; 0; 0; 0];    % initial state: arm slightly offset, pendulum hanging
x_ref  = [pi; 0; pi; 0];    % target: arm at pi, pendulum upright

% Cost weights 
W_x = diag([10, 1, 100, 1]);
W_u = 0.01;
W_T = 10*W_x;   % terminal weight

% new valas
W_x = diag([10, 1, 100, 1]);
W_u = 1e-3;
W_T = 5*W_x;   % terminal weight

%% -----------------------------------------------------------------------
%  CasADi Symbolic Variables
% ------------------------------------------------------------------------
theta1  = SX.sym('theta1');
theta1p = SX.sym('theta1p');
theta2  = SX.sym('theta2');
theta2p = SX.sym('theta2p');
u_sym   = SX.sym('u');

x = vertcat(theta1, theta1p, theta2, theta2p);
u = u_sym;

%% -----------------------------------------------------------------------
%  Smooth Dynamics  (from get_NonLinPendulumModel_Coulomb.m, Coulomb removed)
%
%  Original RHS:
%    rhs1 = (-b1vis - kM^2/RM)*theta1p
%           - hJ2*sin(2*theta2)*theta1p*theta2p
%           + m2*L1*l2*sin(theta2)*theta2p^2
%           + (kM/RM)*kSteller*u
%           - b1coul*tanh(...)   <-- removed, handled by PSS
%
%    rhs2 = -b2vis*theta2p
%           + 0.5*hJ2*sin(2*theta2)*theta1p^2
%           - l2*m2*sin(theta2)*grav
%           - b2coul*tanh(...)   <-- removed, handled by PSS
%
%    M * [ddtheta1; ddtheta2] = [rhs1_smooth; rhs2_smooth]
%    M = [hJ0 + hJ2*sin^2(theta2),  m2*L1*l2*cos(theta2)]
%        [m2*L1*l2*cos(theta2),      hJ2                  ]
% ------------------------------------------------------------------------
s2     = sin(theta2);
c2     = cos(theta2);
s2_2   = sin(2*theta2);     % = 2*sin(theta2)*cos(theta2)
m2L1l2 = m2 * L1 * l2;

MassMa = [hJ0 + hJ2*s2^2,  m2L1l2*c2;
          m2L1l2*c2,        hJ2       ];

rhs1_s = (-b1vis - kM^2/RM) * theta1p ...
         - hJ2*s2_2 * theta1p * theta2p ...
         + m2L1l2*s2 * theta2p^2 ...
         + (kM/RM)*kSteller * u_sym;

rhs2_s = -b2vis * theta2p ...
         + 0.5*hJ2*s2_2 * theta1p^2 ...
         - l2*m2*s2 * grav;

ddq_s  = inv(MassMa)*vertcat(rhs1_s, rhs2_s);

% Smooth base: xdot = [theta1p; ddtheta1; theta2p; ddtheta2]
f_base = vertcat(theta1p, ddq_s(1), theta2p, ddq_s(2));

%% -----------------------------------------------------------------------
%  Coulomb PSS — additive acceleration increments via M^{-1}
%
%  Coulomb forces enter the RHS BEFORE M^{-1}:
%    rhs1 += -b1coul * sign(theta1p)   (arm joint)
%    rhs2 += -b2coul * sign(theta2p)   (pendulum joint)
%
%  State-space increment from arm Coulomb (theta1p > 0):
%    d_rhs = [-b1coul; 0]   =>   d_ddq = M \ [-b1coul; 0]
%    => df_state = [0; d_ddq(1); 0; d_ddq(2)]   (symbolic in theta2!)
%
%  Because M is off-diagonal, arm Coulomb also perturbs ddtheta2,
%  and pendulum Coulomb also perturbs ddtheta1.
%
%  This is the key difference from the acrobot, where friction enters
%  as a direct generalised force (no M^{-1} needed for the columns).
% ------------------------------------------------------------------------

% Arm Coulomb increments (sign depends on theta1p; magnitude flips)
ddq_c1_pos = inv(MassMa)*[-b1coul; 0];   % theta1p > 0: sign(theta1p)=+1
ddq_c1_neg = inv(MassMa)*[ b1coul; 0];   % theta1p < 0: sign(theta1p)=-1

% Pendulum Coulomb increments
ddq_c2_pos  =inv(MassMa)*[0; -b2coul];   % theta2p > 0
ddq_c2_neg = inv(MassMa)*[0;  b2coul];   % theta2p < 0

% Full state-space increments [0; d_ddtheta1; 0; d_ddtheta2]
df_c1_pos = vertcat(0, ddq_c1_pos(1), 0, ddq_c1_pos(2));
df_c1_neg = vertcat(0, ddq_c1_neg(1), 0, ddq_c1_neg(2));
df_c2_pos = vertcat(0, ddq_c2_pos(1), 0, ddq_c2_pos(2));
df_c2_neg = vertcat(0, ddq_c2_neg(1), 0, ddq_c2_neg(2));

% Simplifed friction model
% df_c1_pos = vertcat(0, ddq_c1_pos(1), 0, 0);
% df_c1_neg = vertcat(0, ddq_c1_neg(1), 0, 0);
% df_c2_pos = vertcat(0, 0, 0, ddq_c2_pos(2));
% df_c2_neg = vertcat(0, 0, 0, ddq_c2_neg(2));

% PSS mode matrices — ONLY the Coulomb correction (f_base goes into model.f_0)
F1 = [df_c1_pos, df_c1_neg];   % [4x2]: col1 for theta1p>0, col2 for theta1p<0
F2 = [df_c2_pos, df_c2_neg];   % [4x2]: col1 for theta2p>0, col2 for theta2p<0

S1 = [1; -1];   % sign matrix: row1 active for c>0, row2 for c<0
S2 = [1; -1];

c1_sw = theta1p;   % switching surface 1
c2_sw = theta2p;   % switching surface 2

%% -----------------------------------------------------------------------
%  Cost Function
%
% ------------------------------------------------------------------------
% err_vec = [pi*(1+cos(x(1)/2)); x(2); pi*(1+cos(x(3)/2)); x(4)] - x_ref;
err_vec = x-x_ref;
f_q   = err_vec' * W_x * err_vec  +  W_u * u_sym^2;
f_q_T = err_vec' * W_T * err_vec;

%% -----------------------------------------------------------------------
%  nosnoc Problem Options
% ------------------------------------------------------------------------
problem_options = nosnoc.Options();
problem_options.rk_scheme              = RKSchemes.RADAU_IIA;
problem_options.n_s                    = n_s;
problem_options.cross_comp_mode        = "FE_FE";
problem_options.N_stages               = N_stages;
problem_options.N_finite_elements      = N_FE;
problem_options.T                      = T_pred;
problem_options.use_fesd               = use_fesd;
problem_options.gamma_h                = gamma_h;
problem_options.euler_cost_integration = 0;

%% -----------------------------------------------------------------------
%  Homotopy Solver Options
% ------------------------------------------------------------------------
homotopy_options = nosnoc.reg_homotopy.Options();
homotopy_options.complementarity_tol             = tol;
homotopy_options.N_homotopy                      = N_homotopy;
homotopy_options.print_level                     = 3;
homotopy_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
homotopy_options.homotopy_steering_strategy      = "DIRECT";
homotopy_options.lift_complementarities          = 0;
homotopy_options.opts_casadi_nlp.ipopt.max_iter  = 1000;

%% -----------------------------------------------------------------------
%  CCOpt Solver Options
% ------------------------------------------------------------------------
ccopt_options = nosnoc.ccopt.Options(); % ccopt options
% ccopt options
ccopt_options.opts_madnlp.linear_solver = 'Ma27Solver';
ccopt_options.opts_madnlp.tol = tol;
ccopt_options.opts_madnlp.barrier.TYPE = 'MonotoneUpdate';
ccopt_options.opts_madnlp.disable_garbage_collector = true;
% ccopt_options.opts_madnlp.barrier.mu_min = 1e-9;
% ccopt_options.opts_ccopt.relaxation_update.TYPE = 'RolloffRelaxationUpdate';
% ccopt_options.opts_ccopt.relaxation_update.rolloff_slope = 2.5;
% ccopt_options.opts_ccopt.relaxation_update.rolloff_point = 1e-12;
% ccopt_options.opts_ccopt.sigma_min = 1e-8;
%ccopt_options.opts_ccopt.relaxation_update.TYPE = 'RelaxLBUpdate';
%ccopt_options.opts_ccopt.relaxation_update.relax_threshold = 1e-5;
ccopt_options.opts_ccopt.q_regularization = 'critical_rho';
ccopt_options.opts_ccopt.critical_rho_factor = 0.9999;
%ccopt_options.opts_madnlp.barrier.mu_init = 1.0;
%ccopt_options.opts_madnlp.barrier.TYPE = 'QualityFunctionUpdate';

%% -----------------------------------------------------------------------
%  MPC Options
% ------------------------------------------------------------------------
mpc_options = nosnoc.mpc.Options();
mpc_options.solve_advanced_problem          = 0;
mpc_options.advanced_n_qpecs               = 1;
mpc_options.discard_constraints_in_hessian = true;
mpc_options.sqpec_hessian_convexification  = "PROJECT";
mpc_options.use_probing_qp                 = 0;
mpc_options.objective_ratio                = 0.99;

%% -----------------------------------------------------------------------
%  nosnoc Model
% ------------------------------------------------------------------------
model = nosnoc.model.Pss();

model.x   = x;
model.x0  = x0_val;
model.lbx = [-inf; -inf; -inf; -inf];
model.ubx = [ inf;  inf;  inf;  inf];

model.u   = u;
model.lbu = -U_max;
model.ubu =  U_max;

model.f_0 = f_base;    % smooth dynamics (separated to reduce nonlinearity)
model.S   = {S1, S2};
model.c   = {c1_sw, c2_sw};
model.F   = {F1, F2};

% model.f_q   = f_q;
% model.f_q_T = f_q_T;
model.lsq_x = {x, x_ref, W_x};
model.lsq_u = {u, 0, W_u};
model.lsq_T = {x, [x_ref], W_T};

%% -----------------------------------------------------------------------
%  Optional Plant Integrator (model_plant_mismatch = 1)
% ------------------------------------------------------------------------
if model_plant_mismatch
    sim_model       = nosnoc.model.Pss();
    sim_model.x     = x;         sim_model.x0  = x0_val;
    sim_model.lbx   = -inf(4,1); sim_model.ubx =  inf(4,1);
    sim_model.u     = u;         sim_model.lbu = -inf;  sim_model.ubu = inf;
    sim_model.f_0   = f_base;
    sim_model.S     = {S1, S2};  sim_model.c = {c1_sw, c2_sw};
    sim_model.F     = {F1, F2};

    sim_problem_options                    = nosnoc.Options();
    sim_problem_options.rk_scheme          = RKSchemes.RADAU_IIA;
    sim_problem_options.n_s                = n_s+1;
    sim_problem_options.N_finite_elements  = 2;
    sim_problem_options.T_sim              = DT;
    sim_problem_options.N_sim              = 5;
    sim_problem_options.print_level        = 0;

    integrator_options             = nosnoc.integrator.Options();
    integrator_options.print_level = 1;
    sim_solver_opts                = integrator_options.fesd_solver_opts;
    sim_solver_opts.print_level                = 0;
    sim_solver_opts.homotopy_steering_strategy = HomotopySteeringStrategy.ELL_INF;
    sim_solver_opts.complementarity_tol        = 1e-10;

    integrator = nosnoc.Integrator(sim_model, sim_problem_options, integrator_options);
end

%% -----------------------------------------------------------------------
%  Set up MPC Solver (real-time or fully converged, ccopt or not)
% ------------------------------------------------------------------------
if use_rtmpc
    if use_ccopt
        mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, ccopt_options, ccopt_options);
    else
        mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, homotopy_options, homotopy_options);
    end
    % mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, homotopy_options, mpecopt_options);
else
    if use_ccopt
        mpc = nosnoc.mpc.FullMpcc(model, mpc_options, problem_options, ccopt_options);
    else
        mpc = nosnoc.mpc.FullMpcc(model, mpc_options, problem_options, homotopy_options);
    end
end

%% -----------------------------------------------------------------------
%  Live Plot Setup
% ------------------------------------------------------------------------
if plot_intermediate
    fig_live = figure('Name','Furuta MPC – Live','Color','w','Position',[50 50 950 700]);

    ax_q = subplot(3,1,1); hold on; grid on;
    xlim([0 (N_steps + N_stages)*DT]); ylim([-pi*1.15, pi*1.15]);
    ylabel('$\theta$ [rad]'); xlabel('$t$ [s]');
    title('Arm ($\theta_1$, blue) and Pendulum ($\theta_2$, red) Angles');

    ax_v = subplot(3,1,2); hold on; grid on;
    xlim([0 (N_steps + N_stages)*DT]); ylim([-15, 15]);
    ylabel('$\dot\theta$ [rad/s]'); xlabel('$t$ [s]');
    title('Angular Velocities');

    ax_u = subplot(3,1,3); hold on; grid on;
    xlim([0 (N_steps + N_stages)*DT]); ylim([-U_max*1.15, U_max*1.15]);
    ylabel('$u$ [V]'); xlabel('$t$ [s]');
    title('Control Input');
    yline(ax_u,  U_max, 'r--', 'LineWidth', 1.0);
    yline(ax_u, -U_max, 'r--', 'LineWidth', 1.0);

    sgtitle('Furuta Pendulum MPC', 'FontSize',13,'FontWeight','bold');
end

%% -----------------------------------------------------------------------
%  MPC Loop
% ------------------------------------------------------------------------
x_hist = x0_val;   % accumulate closed-loop states
u_hist = [];
t_hist = 0;
x0     = x0_val;

preparation_times = [];
feedback_times    = [];

fprintf('Starting MPC loop (%d steps, T_sim=%.1f s)...\n', N_steps, T_sim);

for step = 1:N_steps

    % --- Preparation ---
    stats_prep = mpc.do_preparation();
    preparation_times(end+1) = stats_prep.preparation_time;

    % --- Feedback ---
    [u_i, stats_fb] = mpc.get_feedback(x0);
    feedback_times(end+1) = stats_fb.feedback_time;

    fprintf('Step %3d/%d | u=%7.3f V | theta2=%.3f rad | prep=%.3fs fb=%.3fs\n', ...
        step, N_steps, u_i, x0(3), stats_prep.preparation_time, stats_fb.feedback_time);

    % --- Advance plant ---
    if model_plant_mismatch
        integrator.set_x0(x0);
        [~, x_int] = integrator.simulate("u", repmat(u_i,[1,5]), "x0", x0);
        x0 = x_int(:,end);
    else
        x0 = mpc.get_predicted_state();
    end

    % --- Store ---
    x_hist = [x_hist, x0];
    u_hist = [u_hist, u_i];
    t_hist = [t_hist, t_hist(end) + DT];

    % --- Live plot ---
    if plot_intermediate
        x_pred      = mpc.get('x');
        u_pred      = mpc.get('u');
        t_grid_pred = mpc.get_time_grid();
        t_grid_u    = mpc.get_control_grid();

        cla(ax_q); hold(ax_q,'on'); grid(ax_q,'on');
        plot(ax_q, t_hist, x_hist(1,:), 'b-',  'LineWidth',1.5);
        plot(ax_q, t_hist, x_hist(3,:), 'r-',  'LineWidth',1.5);
        plot(ax_q, t_hist(end)+t_grid_pred, x_pred(1,:), 'b--','LineWidth',0.9);
        plot(ax_q, t_hist(end)+t_grid_pred, x_pred(3,:), 'r--','LineWidth',0.9);
        yline(ax_q, x_ref(1), 'b:', 'LineWidth',0.8);
        yline(ax_q, x_ref(3), 'r:', 'LineWidth',0.8);
        legend(ax_q,{'$\theta_1$','$\theta_2$','$\theta_1$ pred','$\theta_2$ pred'},...
               'Location','best','FontSize',8);
        ylabel(ax_q,'$\theta$ [rad]'); xlabel(ax_q,'$t$ [s]');

        cla(ax_v); hold(ax_v,'on'); grid(ax_v,'on');
        plot(ax_v, t_hist, x_hist(2,:), 'b-',  'LineWidth',1.5);
        plot(ax_v, t_hist, x_hist(4,:), 'r-',  'LineWidth',1.5);
        plot(ax_v, t_hist(end)+t_grid_pred, x_pred(2,:), 'b--','LineWidth',0.9);
        plot(ax_v, t_hist(end)+t_grid_pred, x_pred(4,:), 'r--','LineWidth',0.9);
        yline(ax_v, x_ref(2), 'b:', 'LineWidth',0.8);
        yline(ax_v, x_ref(4), 'r:', 'LineWidth',0.8);
        ylabel(ax_v,'$\dot\theta$ [rad/s]'); xlabel(ax_v,'$t$ [s]');

        cla(ax_u); hold(ax_u,'on'); grid(ax_u,'on');
        if ~isempty(u_hist)
            stairs(ax_u, t_hist(1:end-1), u_hist, 'k-', 'LineWidth',1.5);
        end
        stairs(ax_u, t_hist(end)+t_grid_u, [u_pred(:)', u_pred(end)], 'k--','LineWidth',0.9);
        yline(ax_u,  U_max,'r--','LineWidth',1.0);
        yline(ax_u, -U_max,'r--','LineWidth',1.0);
        ylabel(ax_u,'$u$ [V]'); xlabel(ax_u,'$t$ [s]');

        drawnow;
    end
end

fprintf('\nMPC loop complete.\n');
fprintf('Final: theta1=%.3f rad (ref pi=%.3f) | theta2=%.3f rad (ref pi=%.3f)\n',...
    x_hist(1,end), pi, x_hist(3,end), pi);

%% -----------------------------------------------------------------------
%  Resulting closed-Loop trajectories
% ------------------------------------------------------------------------
figure('Name','Furuta MPC – Closed-Loop Result','Color','w','Position',[100 100 850 700]);
clrs = [0.18 0.50 0.90; 0.15 0.72 0.48; 0.92 0.38 0.18; 0.70 0.28 0.82];
state_labels = {'$\theta_1$ [rad]','$\dot\theta_1$ [rad/s]','$\theta_2$ [rad]','$\dot\theta_2$ [rad/s]'};
for i = 1:4
    subplot(5,1,i);
    plot(t_hist, x_hist(i,:),'LineWidth',2.0,'Color',clrs(i,:)); hold on;
    yline(x_ref(i),'k--','LineWidth',1.0);
    ylabel(state_labels{i}); 
    xlabel('$t$ [s]');
    grid on; box on;
    xlim([0 T_sim])
    set(gca,'FontSize',12)
end
subplot(5,1,5);
stairs(t_hist(1:end), [u_hist,nan],'LineWidth',2.0,'Color',[0.20 0.20 0.75]); hold on;
yline( U_max,'r--','LineWidth',1.0); yline(-U_max,'r--','LineWidth',1.0);
ylabel('$u$ [V]'); 
xlabel('$t$ [s]');
xlim([0 T_sim])
ylim([-1.1*U_max 1.1*U_max])
set(gca,'FontSize',12)
grid on; box on;
sgtitle('Furuta Pendulum MPC – Closed Loop','FontWeight','bold','FontSize',13);

%% -----------------------------------------------------------------------
%  Timing figure
% ------------------------------------------------------------------------
figure('Name','MPC Timing','Color','w','Position',[100 50 700 300]);
stairs(preparation_times,'DisplayName','Preparation [s]','LineWidth',2); hold on;
stairs(feedback_times,   'DisplayName','Feedback [s]',   'LineWidth',2);
legend show; grid on;
ylabel('Time [s]'); xlabel('MPC step');
title('MPC Step Computation Times');

%% -----------------------------------------------------------------------
%  Animation: 3D pendulum 
% ------------------------------------------------------------------------
fprintf('Generating animation...\n');

th1_cl  = x_hist(1,:);
th1p_cl = x_hist(2,:);
th2_cl  = x_hist(3,:);
th2p_cl = x_hist(4,:);
N_cl    = length(t_hist);

% Figure sized for a clean 16:9 video frame
fig_anim = figure('Name','Furuta Pendulum – Animation',...
                  'Color',[0.05 0.05 0.10],...
                  'Position',[200 100 1280 720]);

ax3 = axes('Parent', fig_anim);
set(ax3,'Color',[0.05 0.05 0.10],...
    'XColor',[0.65 0.65 0.65],'YColor',[0.65 0.65 0.65],'ZColor',[0.65 0.65 0.65],...
    'GridColor',[0.22 0.22 0.22],'GridAlpha',0.6,'FontSize',11,...
    'Position',[0.05 0.05 0.90 0.88]);
hold(ax3,'on'); grid(ax3,'on'); axis(ax3,'equal');

R = (L1 + l2) * 1.4;
xlim(ax3,[-R R]); ylim(ax3,[-R R]); zlim(ax3,[-l2*1.6, l2*1.6]);
xlabel(ax3,'$x$ [m]','Color',[0.8 0.8 0.8],'FontSize',12);
ylabel(ax3,'$y$ [m]','Color',[0.8 0.8 0.8],'FontSize',12);
zlabel(ax3,'$z$ [m]','Color',[0.8 0.8 0.8],'FontSize',12);
title(ax3,'Furuta Pendulum – MPC Swing-Up', 'Color','w','FontSize',13,'FontWeight','bold');
view(ax3, 38, 22);

% Motor axis pillar
plot3(ax3,[0 0],[0 0],[-l2*1.4 l2*1.4],'Color',[0.50 0.50 0.50],'LineWidth',2);
scatter3(ax3,0,0,0,100,'w','filled');

% Ghost: target configuration (theta1=pi, theta2=pi)
th1g = x_ref(1); th2g = x_ref(3);
Axg  = L1*cos(th1g); Ayg = L1*sin(th1g);
Bxg  = Axg + l2*sin(th2g)*(-sin(th1g));
Byg  = Ayg + l2*sin(th2g)*( cos(th1g));
Bzg  = -l2*cos(th2g);
plot3(ax3,[0 Axg],[0 Ayg],[0 0],'--','Color',[0.25 0.55 0.25],'LineWidth',1.5);
plot3(ax3,[Axg Bxg],[Ayg Byg],[0 Bzg],'--','Color',[0.25 0.80 0.25],'LineWidth',1.5);
scatter3(ax3, Bxg, Byg, Bzg, 120,'g','filled');
text(ax3, Bxg+0.006, Byg, Bzg+0.010,'target',...
     'Color',[0.4 0.9 0.4],'FontSize',10,'FontWeight','bold');

% Animated handles
h_arm   = plot3(ax3,[0 0],[0 0],[0 0],'-','Color',[0.35 0.65 1.00],'LineWidth',5);
h_jnt   = scatter3(ax3,0,0,0,   80,[0.35 0.65 1.00],'filled');
h_rod   = plot3(ax3,[0 0],[0 0],[0 0],'-','Color',[1.00 0.55 0.15],'LineWidth',3.5);
h_bob   = scatter3(ax3,0,0,0,  160,[1.00 0.55 0.15],'filled');
h_trail = plot3(ax3,[],[],[],':','Color',[1.00 0.55 0.15],'LineWidth',1.2);

% Time and state labels (top-left of axes)
h_tlbl  = text(ax3,-R*0.93,-R*0.93, l2*1.45, 't = 0.00 s',...
               'Color','w','FontSize',12,'FontWeight','bold');
h_slbl  = text(ax3,-R*0.93,-R*0.93, l2*1.20,...
               '\theta_2 = 0.000 rad',...
               'Color',[1.0 0.7 0.3],'FontSize',10,'Interpreter','tex');

% Open video writer
if save_video
    vid = VideoWriter('furuta_mpc.mp4','MPEG-4');
    vid.FrameRate = round(1/DT);  % real-time
    open(vid);
end

% Animation loop
trail_x = []; trail_y = []; trail_z = [];

for k = 1:N_cl
    t1 = th1_cl(k);
    t2 = th2_cl(k);

    Ax = L1 * cos(t1);
    Ay = L1 * sin(t1);

    Bx = Ax + l2 * sin(t2) * (-sin(t1));
    By = Ay + l2 * sin(t2) * ( cos(t1));
    Bz =    - l2 * cos(t2);

    set(h_arm,  'XData',[0, Ax],  'YData',[0, Ay],  'ZData',[0,  0]);
    set(h_jnt,  'XData', Ax,      'YData', Ay,       'ZData', 0);
    set(h_rod,  'XData',[Ax, Bx], 'YData',[Ay, By], 'ZData',[0, Bz]);
    set(h_bob,  'XData', Bx,      'YData', By,       'ZData', Bz);

    trail_x(end+1) = Bx; 
    trail_y(end+1) = By; 
    trail_z(end+1) = Bz; 
    set(h_trail,'XData',trail_x,'YData',trail_y,'ZData',trail_z);

    set(h_tlbl,'String', sprintf('t = %.2f s', t_hist(k)));
    set(h_slbl,'String', sprintf('\\theta_2 = %.3f rad  (target \\pi = %.3f)', t2, pi));

    drawnow;

    if save_video
        frame = getframe(fig_anim);
        writeVideo(vid, frame);
    end

    pause(0.01);
end

if save_video
    close(vid);
    fprintf('Video saved: furuta_mpc.mp4\n');
end
