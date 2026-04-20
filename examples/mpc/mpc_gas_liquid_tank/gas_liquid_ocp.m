clear; clc; close all;
import casadi.*
%% Description
% Mpc for An ideal gas-liquid tank system with a pressure
% control valve
% Example from: 
% Globally Convergent Method for Optimal Control of Hybrid Dynamical Systems by  Saif R. Kazi, Mandar Thombre,Lorenz Biegler 
% 12th IFAC International Symposium on Advanced Control of Chemical Processes July 14-17, 2024. Toronto, Canada

%% plot and video settings
latexify_plot()
nice_plot_colors;

%%
filename_pdf = 'liquid_gas_mpc';
linewidth = 1.5;
use_rtmpc = 1;

%% Model settings
N_stages = 100; % number of control stages;
DT = 0.25;
N_FE = 3;
use_fesd = 1;
gamma_h = 1.0;

%% solver options
tol = 1e-6;
N_homotopy = 6;
sigma_0 = 1;
fullmpcc_fast_sigma_0 = 1e-3; % sigma for mpec warmstart in mpc

%% nosnoc setup
problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
solver_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosonc problems.

%% Model variables and parameters
% specify initial and desired state

% CasADi symbolic variables
n_x = 2;
n_u = 1;
M_G = SX.sym('M_G'); % rate of vapor holdup
M_L = SX.sym('M_L'); % rate of liquid holdup

x = [M_G;M_L];
u = SX.sym('u'); % control
% Parameters
F_L = 2.5; % mol/sec
F_G = 0.1; % mol/sec
V = 10;
Vs = 5; % liters
T = 300;% K;
P_out = 1; % atm
rho_L = 50; 
k_L = 1;
k_G = 1;
% inital vals;
L0 = 0.25;
G0 = 0.25;
P0 = 35;
M_L0 = 260;
M_G0 = 6.83;
% R = 8.314; % J/mol K
R = 0.0821; %L atm/(mol K)

%% Model

c = M_L/rho_L-Vs;

% ideal gas equation
% M_G*R*T/P+ M_L/rho_L = V;
% M_G*R*T/P= V-M_L/rho_L;
% P/M_G*R*T= 1/V-M_L/rho_L;
P= M_G*R*T/(V-M_L/rho_L);

P_fun = Function('P_fun', {x}, {P});
c_fun = Function('P_fun', {x}, {c});

L = k_L*u*(P-P_out);
G = k_G*u*(P-P_out);

% c >0  liquid mode
f_1 = [F_G;F_L-L];

% c < 0  liquid mode
f_2 = [F_G-G; F_L];

ubx = [inf; inf];
lbx = [-inf; -inf];
u_max = 0.15;
u_ref = 0.1;

% stage cost
f_q = 100*(u-u_ref)^2;
% terminal cost
f_q_T  = M_L;
F = [f_1, f_2];
x0 = [M_G0;M_L0];
% Populate nosnoc model
model = nosnoc.model.Pss();
model.c = c;         % switching function c: cart velocity
model.S = [1; -1];   % sign matrix S % f_1 for c>0, f_2 for c<0
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
problem_options.T = N_stages*DT;  % Time horizon
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 2;
problem_options.dcs_mode = 'Stewart';
problem_options.N_stages = N_stages; % number of control intervals
problem_options.N_finite_elements = N_FE; % number of finite element on every control interval
problem_options.cross_comp_mode = "FE_FE";
problem_options.use_fesd = use_fesd;
problem_options.gamma_h = gamma_h;

%% MPEC solver options
solver_options.complementarity_tol = tol; % Value to drive the complementarity residual to.
solver_options.N_homotopy = N_homotopy; % Maximum number of homotopy iterations.
solver_options.print_level = 3;
solver_options.sigma_0 = 1;
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; %
solver_options.homotopy_steering_strategy = "ELL_INF";

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
stairs(t_u,[u,u(end)],'LineWidth',linewidth)
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

