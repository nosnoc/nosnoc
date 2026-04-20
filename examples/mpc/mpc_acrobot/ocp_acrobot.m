%% Acrobot with joint friction.
% solve single optimal control problem for swing up.
%%
clear all; clc;
close all;
%%
import casadi.*
import nosnoc.*
latexify_plot()
%% Basic timing settings
N_FE = 3;
use_fesd = 1;
gamma_h = 1;
n_s = 2;
dt = 0.1;  % Time step
T = 4;    % Total simulation time
N_stages = T/dt;   % Number of steps

%% Define model parameters
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
% F_c = [-mu1*tanh(10*q1_dot); -mu2*tanh(10*q2_dot)];  % Coulomb friction

% Compute accelerations
q_ddot = M \ (u - C - G - F_v);

% Define system dynamics
f = [q1_dot; q2_dot; q_ddot];

% Create CasADi function

% Define bounds for states and controls
lbx = [-pi; -pi; -inf; -inf]; % Lower bound on states
ubx = [pi; pi; inf; inf];    % Upper bound on states
lbu = [-10; -10];             % Lower bound on controls
ubu = [10; 10];               % Upper bound on controls

% Define reference trajectories
x_ref = [pi; 0; 0; 0];  % Upright position
u_ref = [0; 0];        % Zero control
x0 = [0;0;0;0];

% Define weight matrices for cost function
Q_cost = diag([10, 10, 1, 1]);   % State cost weights
R_cost = diag([0.1, 0.1]);       % Control cost weights
P_cost = 10*Q_cost;


%% all options
problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
solver_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosonc problems.

%%
% Choosing the Runge - Kutta Method and number of stages
problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
problem_options.n_s = n_s; % Number of stage points in the RK method (determines accuracy)
problem_options.cross_comp_mode = "FE_FE";
% Time-settings  - Solve a
problem_options.N_stages = N_stages; % Number of control intervals over the time horizon.
problem_options.N_finite_elements = N_FE; % Number of finite element (integration steps) on every control interval (optionally a vector might be passed).
problem_options.T = T;    % Time horizon (for a time-optimal problem, the actual horizon length is obtained by solving the problem).
problem_options.use_fesd = use_fesd;
problem_options.gamma_h = gamma_h;
%% nosnoc models
model = nosnoc.model.Pss(); %
model.x = x;
model.x0 = x0;
model.lbx = lbx; % lower bounds on states
model.ubx = ubx; % upper bounds on states
% define control vectors
model.u = u;
model.lbu = lbu;
model.ubu = ubu;

% Equation of motion
q_ddot = M \ (u - C - G - F_v);

% Define dynamics
f_base = [q1_dot; q2_dot; q_ddot];

f_11 = f_base + [0;0;-mu1;0];
f_12 = f_base + [0;0;mu1;0];

f_21 = [0;0;0;-mu2];
f_22 = [0;0;0;mu2];

c1 = q1_dot;
c2 = q2_dot;
% sign matrix for the modes
S1 = [1;-1];
S2 = [1;-1];
F1 = [f_11 f_12];
F2 = [f_21 f_22];

model.S = {S1,S2};
model.c = {c1,c2};
model.F = {F1 F2};

model.f_q = (x-x_ref)'*Q_cost*(x-x_ref) + (u-u_ref)'*R_cost*(u-u_ref); % Add stage cost
model.f_q_T = (x-x_ref)'*P_cost*(x-x_ref); % Add terminal quadratic cost

%% Setup solver an ocp 
solver_options.complementarity_tol = 1e-6;
solver_options.N_homotopy = 10;
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'mumps'; % Using an HSL solver is a significant speedup over the default 'mumps' but requires installation
solver_options.homotopy_steering_strategy ="DIRECT";
solver_options.lift_complementarities = 0;

%% Create and solve OCP
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();

%% plots
% unfold structure to workspace of this script
x = ocp_solver.get('x');
u = ocp_solver.get("u");
t = ocp_solver.get_time_grid();
t_u = ocp_solver.get_control_grid();

%% Plot
figure
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
stairs(t_u,[u,u(:,end)]','LineWidth',1.5)
hold on
% yline(u_max,'k--','LineWidth',1.5)
% yline(-u_max,'k--','LineWidth',1.5)
xlabel("$t$")
ylabel("$u$")
% ylim([-1.1*u_max 1.1*u_max])


