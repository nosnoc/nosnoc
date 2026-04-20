% Solve an optimal control problem with an osc. with variable stiffnes.

clear; clc; close all;
import casadi.*
import nosnoc.*
%% Parameters
N_stages = 100;
m = 1.0;       % Mass (kg)
k1 = 0.5;      % Spring constant in negative region (N/m)
k2 = 3.0;      % Spring constant in positive region (N/m)
c = 0.7;       % Damping coefficient (Ns/m)

Q = diag([1.0, 0.5]);  % State tracking weight
R = 0.1;               % Control input weight
P = diag([5.0, 2.5]);  % Terminal cost weight

x1_max = 3.0;
x2_max = 4.0;
% Input constraints
x0 = [-1.0; 0.0];
x_ref = [1.5; 0.0];    % Reference state [position; velocity]
% x_ref = [0; 0.0];    % Reference state [position; velocity]

u_ss1 = -k1/m*x_ref(1) - c/m*x_ref(2);
u_ss2 = -k2/m*x_ref(1) - c/m*x_ref(2);
u_max = 5; % set point reachable
% u_max = 2; % set point unreachable
u_ref = -1*u_ss1;
u_ref = -u_ss2;


%% all options
problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
solver_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosonc problems.

%% problem options and model
% Choosing the Runge - Kutta Method and number of stages
problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
problem_options.n_s = 3; % Number of stage points in the RK method (determines accuracy)
problem_options.cross_comp_mode = "FE_FE";

% Time-settings  - Solve a
problem_options.N_stages = N_stages; % Number of control intervals over the time horizon.
problem_options.N_finite_elements = 2; % Number of finite element (integration steps) on every control interval (optionally a vector might be passed).
problem_options.T = 5;    % Time horizon (for a time-optimal problem, the actual horizon length is obtained by solving the problem).

model = nosnoc.model.Pss(); % Initialize a nosnoc model, (depending on the underlying nonsmooth model, several options exist)


%% model
% define differential states and populate the model.
x = SX.sym('x',2);
u = SX.sym('u');  % CasADi symbolic variables for controls
model.x = x;
model.x0 = x0;
model.lbx = [-x1_max;-x2_max]; % lower bounds on states
model.ubx = [x1_max;x2_max]; % upper bounds on states
% define control vectors
model.u = u;
model.lbu = -u_max;
model.ubu = u_max;

rhs_neg = [x(2); 
           -k1/m*x(1) - c/m*x(2) + 1/m*u];
rhs_pos = [x(2); 
           -k2/m*x(1) - c/m*x(2) + 1/m*u];

model.c = x(1); 
model.S = [-1;1]; % Region R_1 is defined by c<0 (hence the -1), region R_2 by c>0 (hence the +1) in the sign matrix S.
model.F = [rhs_neg rhs_pos]; % The columns of this matrix store the vector fields of every region.
model.f_q = (x-x_ref)'*Q*(x-x_ref) + R*(u-u_ref)^2; % Add stage cost
model.f_q_T = (x-x_ref)'*P*(x-x_ref); % Add terminal quadratic cost

% Setup solver an mpc options
homotopy_options.complementarity_tol = 1e-6; 
homotopy_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; % Using an HSL solver is a significant speedup over the default 'mumps' but requires installation
homotopy_options.lift_complementarities = 0;

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
stairs(t_u,[u,u(end)],'LineWidth',1.5)
hold on
yline(u_max,'k--','LineWidth',1.5)
yline(-u_max,'k--','LineWidth',1.5)
xlabel("$t$")
ylabel("$u$")
ylim([-1.1*u_max 1.1*u_max])

