%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adapted for nosnoc / original acados model by Jörg Fischer
% https://github.com/Jo-Fischer/acados-STM32-NUCLEO-H745ZI/tree/master
% Description:
%   Swing-up optimal control problem for the Furuta Pendulum using nosnoc.
%   Coulomb friction is treated as a piecewise-smooth system (PSS) - no tanh smoothing.
%
%   States:  x = [theta1, theta1_dot, theta2, theta2_dot]
%     theta1     : actuator arm angle [rad]   (horizontal rotation)
%     theta1_dot : actuator arm angular velocity [rad/s]
%     theta2     : pendulum angle [rad]  (0 = hanging DOWN, pi = upright)
%     theta2_dot : pendulum angular velocity [rad/s]
%
%   Input:   u  motor voltage [V], constrained to [-U_max, U_max]
%
%   Initial state:  x0 = [0.1; 0; 0; 0]  
%     arm at theta1=0.1 rad, pendulum hanging (theta2=0), everything at rest
%
%   Target state (from acados yref = [pi; 0; pi; 0]):
%     theta1 = pi  (arm rotated 180 degrees)
%     theta2 = pi  (pendulum upright)
%
%   Cost (replicates acados NONLINEAR_LS cost exactly):
%     output  y(theta) = pi*(1 + cos(theta/2))
%     yref              = pi
%     Verify:  theta=0  (hanging) -> y = pi*(1+1) = 2*pi   (bad, costs (pi)^2)
%              theta=pi (upright) -> y = pi*(1+0) = pi     (good, cost=0)
%
%   PSS structure — 2 independent switches (additive decomposition):
%     Switch 1: sign(theta1_dot)  ->  arm Coulomb friction  (b1coul)
%     Switch 2: sign(theta2_dot)  ->  pendulum Coulomb friction (b2coul)
%
%   Coulomb forces enter the RHS before solving M*ddq = rhs:
%     rhs1 += -b1coul * sign(theta1p)
%     rhs2 += -b2coul * sign(theta2p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;

import casadi.*
import nosnoc.*
latexify_plot();

%% -----------------------------------------------------------------------
%  Physical Parameters  (from Parameters_FurutaPendulum_Coulomb.m)
% ------------------------------------------------------------------------
grav     = 9.81;
kSteller = 1.0983;   % PWM unit voltage factor
kM       = 0.0236;   % motor torque constant
RM       = 0.5200;   % motor resistance
b1vis    = 6.0700e-04;  % arm viscous friction
b1coul   = 0.0104;      % arm Coulomb friction
hJ0      = 0.0021;      % inertia of arm + hanging pendulum about motor axis
L1       = 0.1215;      % arm length
m2       = 0.0461;      % pendulum mass
l2       = 0.0790;      % distance from pendulum joint to pendulum CoM
b2vis    = 2.6830e-05;  % pendulum viscous friction
b2coul   = 7.5418e-04;  % pendulum Coulomb friction
hJ2      = 3.7054e-04;  % pendulum inertia about pendulum joint

%% -----------------------------------------------------------------------
%  OCP Timing and Weights
% ------------------------------------------------------------------------
T        = 1;      % total horizon [s]
dt       = 0.02;     % control interval [s]
N_stages = T / dt;   % N control intervals
N_FE     = 2;        % finite elements per interval (accuracy)
n_s      = 2;        % Radau IIA stages

U_max    = 10.0;     % [V] input saturation

x0_val   = [0.1; 0; 0; 0];   % initial state 
x_ref = [pi; 0; pi; 0];      % target states

% Cost weights 
W_x = diag([100, 5, 100, 1]);  % state weights
W_u = 0.01;                   % input weight
W_T = 10*W_x;                  % terminal weight 

%% -----------------------------------------------------------------------
%  CasADi Symbolic Variables
% ------------------------------------------------------------------------
theta1  = SX.sym('theta1');    % arm angle
theta1p = SX.sym('theta1p');   % arm angular velocity
theta2  = SX.sym('theta2');    % pendulum angle
theta2p = SX.sym('theta2p');   % pendulum angular velocity
u_sym   = SX.sym('u');         % motor voltage

x = vertcat(theta1, theta1p, theta2, theta2p);
u = u_sym;

%% -----------------------------------------------------------------------
%  We remove the Coulomb terms -> smooth RHS below.
%  The Coulomb terms are handled via PSS modes (see below).
% ------------------------------------------------------------------------
s2     = sin(theta2);
c2     = cos(theta2);
s2_2   = sin(2*theta2);
m2L1l2 = m2 * L1 * l2;

% Mass matrix M(theta2)
MassMa = [hJ0 + hJ2*s2^2,  m2L1l2*c2;
          m2L1l2*c2,        hJ2       ];

% Smooth RHS (viscous damping + gyroscopic/centripetal + gravity + motor)
rhs1_s = (-b1vis - kM^2/RM) * theta1p ...
         - hJ2*s2_2 * theta1p * theta2p ...
         + m2L1l2*s2 * theta2p^2 ...
         + (kM/RM)*kSteller * u_sym;

rhs2_s = -b2vis * theta2p ...
         + 0.5*hJ2*s2_2 * theta1p^2 ...
         - l2*m2*s2 * grav;

% Smooth accelerations
ddq_s  = inv(MassMa)*vertcat(rhs1_s, rhs2_s);

% Smooth base state derivative (no Coulomb)
f_base = vertcat(theta1p, ddq_s(1), theta2p, ddq_s(2));

%% -----------------------------------------------------------------------
%  Coulomb Friction — PSS Acceleration Increments
%
%  Coulomb enters the RHS as:
%    rhs1 += -b1coul * sign(theta1p)    (joint 1, arm)
%    rhs2 += -b2coul * sign(theta2p)    (joint 2, pendulum)
%
%  The resulting acceleration corrections are M \ [coulomb_force]:
%    When theta1p > 0:  rhs1_correction = -b1coul  -> ddq_c1 = M\[-b1coul; 0]
%    When theta1p < 0:  rhs1_correction = +b1coul  -> ddq_c1 = M\[+b1coul; 0]
%  (and analogously for theta2p and b2coul)
%
%  Note: the off-diagonal M coupling means that arm Coulomb also affects
%  ddtheta2, and pendulum Coulomb also affects ddtheta1.
% ------------------------------------------------------------------------

% Arm Coulomb acceleration increments (symbolic, depend on theta2 via M)
ddq_c1_pos = inv(MassMa)*[-b1coul; 0];   % arm Coulomb, theta1p > 0
ddq_c1_neg = inv(MassMa)*[b1coul; 0];   % arm Coulomb, theta1p < 0
% Pendulum Coulomb acceleration increments
ddq_c2_pos = inv(MassMa)*[0; -b2coul];   % pend Coulomb, theta2p > 0
ddq_c2_neg = inv(MassMa)*[0; +b2coul];   % pend Coulomb, theta2p < 0

% Convert to full state-derivative increments [0; dddq1; 0; dddq2]
df_c1_pos = vertcat(0, ddq_c1_pos(1), 0, ddq_c1_pos(2));
df_c1_neg = vertcat(0, ddq_c1_neg(1), 0, ddq_c1_neg(2));
df_c2_pos = vertcat(0, ddq_c2_pos(1), 0, ddq_c2_pos(2));
df_c2_neg = vertcat(0, ddq_c2_neg(1), 0, ddq_c2_neg(2));

%% -----------------------------------------------------------------------
%  PSS Mode Matrices (additive decomposition, same pattern as nosnoc acrobot)
%
%  Subsystem 1 — switch on theta1p (arm Coulomb):
%    F{1}(:,1) = f_base + df_c1_pos   (active when theta1p > 0)
%    F{1}(:,2) = f_base + df_c1_neg   (active when theta1p < 0)
%  Subsystem 2 — switch on theta2p (pendulum Coulomb, additive only):
%    F{2}(:,1) = df_c2_pos            (active when theta2p > 0)
%    F{2}(:,2) = df_c2_neg            (active when theta2p < 0)
%  No f_base here — it is already included in subsystem 1.
%
%  Sign matrices S{k}: S(i)=+1 means column i is active when c_k > 0.
% ------------------------------------------------------------------------
% F1 = [f_base + df_c1_pos,  f_base + df_c1_neg];  % [4 x 2]
% F2 = [df_c2_pos,           df_c2_neg          ];  % [4 x 2]

F1 = [df_c1_pos,  df_c1_neg];  % [4 x 2]
F2 = [df_c2_pos,  df_c2_neg];  % [4 x 2]

S1 = [1; -1];   % col 1 for theta1p>0, col 2 for theta1p<0
S2 = [1; -1];   % col 1 for theta2p>0, col 2 for theta2p<0

c1_sw = theta1p;   % switching surface 1
c2_sw = theta2p;   % switching surface 2

%% -----------------------------------------------------------------------
%  Cost Function  
% ------------------------------------------------------------------------
% err1 = pi * cos(theta1/2);   % arm angle error in output space
% err2 = theta1p;              % arm angular velocity error
% err3 = pi * cos(theta2/2);   % pendulum angle error in output space
% err4 = theta2p;              % pendulum angular velocity error
% err_vec = vertcat(err1, err2, err3, err4);
% err_vec = [pi*(1+cos(x(1)/2)); x(2); pi*(1+cos(x(3)/2)); x(4)]-x_ref;
err_vec = x-x_ref;
f_q   = err_vec' * W_x * err_vec  +  W_u * u_sym^2;
f_q_T = err_vec' * W_T * err_vec;

%% -----------------------------------------------------------------------
%  nosnoc Problem and Solver Options
% ------------------------------------------------------------------------
problem_options = nosnoc.Options();
problem_options.rk_scheme            = RKSchemes.RADAU_IIA;
problem_options.n_s                  = n_s;
problem_options.cross_comp_mode      = "FE_FE";
problem_options.N_stages             = N_stages;
problem_options.N_finite_elements    = N_FE;
problem_options.T                    = T;
problem_options.use_fesd             = 1;
problem_options.gamma_h              = 1;

solver_options = nosnoc.reg_homotopy.Options();
solver_options.complementarity_tol              = 1e-6;
solver_options.N_homotopy                       = 12;
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
solver_options.homotopy_steering_strategy       = "DIRECT";
solver_options.lift_complementarities           = 0;
solver_options.opts_casadi_nlp.ipopt.max_iter = 1e3;

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

model.S = {S1, S2};
model.c = {c1_sw, c2_sw};
model.F = {F1, F2};
model.f_0 = f_base;

model.f_q   = f_q;
model.f_q_T = f_q_T;

%% -----------------------------------------------------------------------
%  Create and Solve OCP
% ------------------------------------------------------------------------
fprintf('Setting up nosnoc OCP solver...\n');
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);

fprintf('Solving swing-up OCP ...\n');
ocp_solver.solve();
fprintf('Done.\n\n');

%% -----------------------------------------------------------------------
%  Extract Solution
% ------------------------------------------------------------------------
x_sol = ocp_solver.get('x');
u_sol = ocp_solver.get('u');
t_x   = ocp_solver.get_time_grid();
t_u   = ocp_solver.get_control_grid();

th1_sol  = x_sol(1,:);
th1p_sol = x_sol(2,:);
th2_sol  = x_sol(3,:);
th2p_sol = x_sol(4,:);

fprintf('--- Solution Summary ---\n');
fprintf('Initial:  theta1=%.3f rad,  theta2=%.3f rad (0=hanging)\n', th1_sol(1),  th2_sol(1));
fprintf('Final:    theta1=%.3f rad  (target: pi = %.3f rad)\n',      th1_sol(end), pi);
fprintf('Final:    theta2=%.3f rad  (target: pi = %.3f rad, upright)\n', th2_sol(end), pi);
fprintf('Final:    theta1_dot=%.4f,  theta2_dot=%.4f rad/s\n', th1p_sol(end), th2p_sol(end));

%% -----------------------------------------------------------------------
%  Figure 1: State and Input Trajectories
% ------------------------------------------------------------------------
figure('Name','Furuta Pendulum – Optimal Trajectories','Color','w','Position',[50 50 850 700]);

state_labels = {'$\theta_1$  [rad]', '$\dot\theta_1$  [rad/s]', '$\theta_2$  [rad]', '$\dot\theta_2$  [rad/s]'};
target_vals  = [pi, 0, pi, 0];
clrs = [0.18 0.50 0.90; 0.15 0.72 0.48; 0.92 0.38 0.18; 0.70 0.28 0.82];

for i = 1:4
    subplot(5,1,i);
    plot(t_x, x_sol(i,:), 'LineWidth', 1.9, 'Color', clrs(i,:)); hold on;
    yline(target_vals(i), 'k--', 'LineWidth', 1.0);
    ylabel(state_labels{i}, 'FontSize', 10);
    xlabel('t [s]', 'FontSize', 10);
    grid on; box on;
end

subplot(5,1,5);
stairs(t_u, [u_sol(:)', u_sol(end)], 'LineWidth', 1.9, 'Color', [0.20 0.20 0.75]);
yline( U_max, 'r--', 'LineWidth', 1.0);
yline(-U_max, 'r--', 'LineWidth', 1.0);
ylabel('u [V]', 'FontSize', 10); xlabel('t [s]', 'FontSize', 9);
grid on; box on;

sgtitle('Furuta Pendulum Swing-Up', ...
        'FontWeight', 'bold', 'FontSize', 13);

%% -----------------------------------------------------------------------
%  Figure 2: Animation
%
%  Geometry (right-hand coordinate system, z positive upward):
%    Motor axis at origin (z-axis).
%    Arm rotates in horizontal plane: arm tip at
%      Ax = L1*cos(theta1),  Ay = L1*sin(theta1),  Az = 0
%
%    Pendulum joint is at the arm tip.
%    Pendulum swings in the vertical plane that contains the arm vector.
%    The out-of-arm-plane horizontal direction (perpendicular to arm):
%      px = -sin(theta1),  py = cos(theta1)
%
%    Rod direction from joint to bob (length l2):
%      horizontal component: l2*sin(theta2)  in the (px,py) direction
%      vertical component:  -l2*cos(theta2)  (negative z when theta2=0)
%
%    Bob position:
%      Bx = Ax + l2*sin(theta2)*(-sin(theta1))
%      By = Ay + l2*sin(theta2)*( cos(theta1))
%      Bz =    - l2*cos(theta2)
%
% ------------------------------------------------------------------------
fprintf('Generating animation...\n');

n_frames = min(300, length(t_x));
idx      = round(linspace(1, length(t_x), n_frames));
t_anim   = t_x(idx);
th1_a    = th1_sol(idx);
th2_a    = th2_sol(idx);
th2p_a   = th2p_sol(idx);

fig_anim = figure('Name','Furuta Pendulum – Animation', ...
                  'Color',[0.05 0.05 0.10], ...
                  'Position',[200 100 1350 740]);

% ---- Left: 3D view ----
ax3 = subplot(1,2,1,'Parent', fig_anim);
set(ax3, 'Color',[0.05 0.05 0.10], ...
    'XColor',[0.65 0.65 0.65], 'YColor',[0.65 0.65 0.65], 'ZColor',[0.65 0.65 0.65], ...
    'GridColor',[0.22 0.22 0.22], 'GridAlpha', 0.6, 'FontSize', 9);
hold(ax3,'on'); grid(ax3,'on'); axis(ax3,'equal');
R = (L1 + l2) * 1.3;
xlim(ax3,[-R R]); ylim(ax3,[-R R]); zlim(ax3,[-l2*1.5, l2*1.5]);
xlabel(ax3,'$x$ [m]','Color',[0.8 0.8 0.8]);
ylabel(ax3,'$y$ [m]','Color',[0.8 0.8 0.8]);
zlabel(ax3,'$z$ [m]','Color',[0.8 0.8 0.8]);
title(ax3,'3D View','Color','w','FontSize',11);
view(ax3, 38, 25);

% Motor axis
plot3(ax3,[0 0],[0 0],[-l2*1.3 l2*1.3],'Color',[0.45 0.45 0.45],'LineWidth',1.5);
scatter3(ax3,0,0,0,90,'w','filled');

% Ghost: target configuration (theta1=pi, theta2=pi)
th1g = pi; th2g = pi;
Axg = L1*cos(th1g); Ayg = L1*sin(th1g);
Bxg = Axg + l2*sin(th2g)*(-sin(th1g));
Byg = Ayg + l2*sin(th2g)*( cos(th1g));
Bzg = -l2*cos(th2g);
plot3(ax3,[0 Axg],[0 Ayg],[0 0],'--','Color',[0.25 0.55 0.25],'LineWidth',1.2);
plot3(ax3,[Axg Bxg],[Ayg Byg],[0 Bzg],'--','Color',[0.25 0.75 0.25],'LineWidth',1.2);
scatter3(ax3,Bxg,Byg,Bzg, 100,'g','filled');
text(ax3,Bxg+0.01,Byg,Bzg+0.01,'target','Color',[0.4 0.9 0.4],'FontSize',8);

% Animated objects
h_arm   = plot3(ax3,[0 0],[0 0],[0 0],'-','Color',[0.35 0.65 1.00],'LineWidth',4.0);
h_jnt   = scatter3(ax3,0,0,0,  70,[0.35 0.65 1.00],'filled');
h_rod   = plot3(ax3,[0 0],[0 0],[0 0],'-','Color',[1.00 0.55 0.15],'LineWidth',3.0);
h_bob   = scatter3(ax3,0,0,0, 150,[1.00 0.55 0.15],'filled');
h_trail = plot3(ax3,NaN,NaN,NaN,':','Color',[1.00 0.55 0.15],'LineWidth',1.0);
h_lbl   = text(ax3,-R*0.92,-R*0.92,l2*1.3,'t = 0.00 s', ...
               'Color','w','FontSize',10,'FontWeight','bold');

% ---- Right: pendulum phase portrait ----
ax_ph = subplot(1,2,2,'Parent', fig_anim);
set(ax_ph,'Color',[0.05 0.05 0.10], ...
    'XColor',[0.65 0.65 0.65],'YColor',[0.65 0.65 0.65], ...
    'GridColor',[0.22 0.22 0.22],'GridAlpha',0.6,'FontSize',9);
hold(ax_ph,'on'); grid(ax_ph,'on');
xlabel(ax_ph,'$\theta_2$  [rad]','Color',[0.8 0.8 0.8],'FontSize',10);
ylabel(ax_ph,'$\dot\theta_2$  [rad/s]','Color',[0.8 0.8 0.8],'FontSize',10);
title(ax_ph,'Pendulum Phase Portrait  $(\theta_2, \dot\theta_2)$','Color','w','FontSize',11);

% Full trajectory (dim background)
plot(ax_ph, th2_sol, th2p_sol,'-','Color',[0.30 0.30 0.30],'LineWidth',0.8);
% Reference lines
xline(ax_ph,  0,'--','Color',[0.40 0.40 0.40],'LineWidth',0.8,'Alpha',0.6);
xline(ax_ph, pi,'--','Color',[0.25 0.65 0.25],'LineWidth',1.0,'Alpha',0.8);
yline(ax_ph,  0,'--','Color',[0.40 0.40 0.40],'LineWidth',0.8,'Alpha',0.6);
% Markers
scatter(ax_ph, th2_sol(1),   th2p_sol(1),   90,'w','o','filled'); % start
scatter(ax_ph, pi,           0,            140,'g','p','filled'); % target

h_ph_trail = plot(ax_ph,NaN,NaN,'-','Color',[1.00 0.55 0.15],'LineWidth',1.8);
h_ph_dot   = scatter(ax_ph, th2_a(1), th2p_a(1), 100,[1.00 0.55 0.15],'filled');

% ---- Animation loop ----
trail_x = []; trail_y = []; trail_z = [];

for k = 1:n_frames
    t1 = th1_a(k);
    t2 = th2_a(k);

    % Arm tip (horizontal plane, z=0)
    Ax = L1 * cos(t1);
    Ay = L1 * sin(t1);

    % Pendulum bob
    %   perpendicular horizontal direction to arm: (-sin(t1), cos(t1))
    Bx = Ax + l2 * sin(t2) * (-sin(t1));
    By = Ay + l2 * sin(t2) * ( cos(t1));
    Bz =    - l2 * cos(t2);     % 0->down(-l2), pi->up(+l2)

    set(h_arm,  'XData',[0, Ax],  'YData',[0, Ay],  'ZData',[0, 0]);
    set(h_jnt,  'XData', Ax,      'YData', Ay,       'ZData', 0);
    set(h_rod,  'XData',[Ax, Bx], 'YData',[Ay, By], 'ZData',[0, Bz]);
    set(h_bob,  'XData', Bx,      'YData', By,       'ZData', Bz);

    trail_x(end+1) = Bx; 
    trail_y(end+1) = By; 
    trail_z(end+1) = Bz; 
    set(h_trail,'XData',trail_x,'YData',trail_y,'ZData',trail_z);

    set(h_ph_trail,'XData',th2_a(1:k),  'YData',th2p_a(1:k));
    set(h_ph_dot,  'XData',th2_a(k),    'YData',th2p_a(k));
    set(h_lbl,'String', sprintf('t = %.2f s', t_anim(k)));

    drawnow limitrate;
    pause(dt/10);
end

fprintf('Animation complete.\n');