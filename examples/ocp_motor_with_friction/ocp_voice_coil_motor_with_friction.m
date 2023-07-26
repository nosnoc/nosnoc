
%% Info
% This is an optimal control example from the paper:
% Optimal control of a voice-coil-motor with Coulombic friction
% by Bahne Christiansen; Helmut Maurer; Oliver Zirn
% Published in: 2008 47th IEEE Conference on Decision and Control
% DOI: 10.1109/CDC.2008.4739025
%% Clear
clc;
clear all; 
close all;
%% Build problem
import casadi.*
problem_options = NosnocProblemOptions();
solver_options = NosnocSolverOptions();
model = NosnocModel();
% Choosing the Runge - Kutta Method and number of stages
problem_options.irk_scheme = IRKSchemes.RADAU_IIA;
problem_options.n_s = 2;
% solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma57';
% MPCC Method
solver_options.N_homotopy = 10;
% Discretization parameters
problem_options.N_stages = 30; % number of control intervals
problem_options.N_finite_elements = 3; % number of finite element on every control intevral (optionally a vector might be passed)
model.T = 0.08;    % Time horizon

%% The Model
% Parameters
m1 = 1.03; % slide mass
m2 = 0.56; % load mass
k = 2.4e3; %  spring constant N/m
c = 0.00; % damping
U_max = 5; % voltage Back-EMF, U = K_s*v_1;
R  = 2; % coil resistance ohm
L = 2e-3; % inductivity, henry
K_F = 12; % force constant N/A ; F_L = K_F*I; % Lorenz force
K_S = 12; % Vs/m (not provided in the paper above)
F_R = 2.1; % guide friction force, N
x0 = [0;0;0;0;0];

%% Symbolic variables
x1 = SX.sym('x1');  % motor mass position
v1 = SX.sym('v1'); % motor mass velocity
x2 = SX.sym('x2');  % load mass position
v2 = SX.sym('v2'); % load mass velocity
I = SX.sym('I'); % electric current
x = [x1;v1;x2;v2;I];
model.x = x;
model.x0 = x0;
% control
U = SX.sym('U'); % the motor voltage
u = [U];
n_u = 1;
model.u = u;
model.lbu = -U_max*ones(n_u,1);
model.ubu = U_max*ones(n_u,1);


%% Dynamics

A = [0  1   0   0   0;...
    -k/m1 -c/m1 k/m1 c/m1 K_F/m1;...
    0   0   0   1   0;...
    k/m2   c/m2   -k/m2   -c/m2 0;...
    0 -K_S/L    0       0   -R/L];
B = [zeros(4,1);1/L];
C1 = [0;-F_R/m1;0;0;0]; % v1 >0
C2 = -C1; %v1<0

% switching dynamics with different friction froces
f_1 = A*x+B*u+C1; % v1>0
f_2 = A*x+B*u+C2; % v1<0

% All modes
F = [f_1, f_2];
% Switching function
c1 = v1;
% Sign matrix (pass a cell when having indepdented subsystems)
model.S = [1;-1];
% The various modes
model.F = F;
% The switching functions
model.c = c1;
% OCP
x_target = [0.01;0;0.01;0;0];
% Stage cost
model.f_q = u^2;
model.g_terminal = x-x_target;


% Inequality constraints
%
% cv = 10; cx = 10;
% model.g_path = [v1-v2;x1-x2];
% model.g_path_ub = [cv;cx];
% model.g_path_lb = -[cv;cx];

%% Solve OCP
mpcc = NosnocMPCC(problem_options, model);
solver = NosnocSolver(mpcc, solver_options);
[results,stats] = solver.solve();

%% plots
% unfold structure to workspace of this script
x1_opt = results.x(1,:);
v1_opt= results.x(2,:);
x2_opt= results.x(3,:);
v2_opt= results.x(4,:);
I_opt= results.x(5,:);

figure
subplot(411)
plot(results.t_grid,x1_opt)
hold on
plot(results.t_grid,x2_opt)
ylabel('$x(t)$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
grid on
legend({'$x_1(t)$','$x_2(t)$'},'Interpreter','latex','Location','best')
subplot(412)
plot(results.t_grid,v1_opt)
hold on
plot(results.t_grid,v2_opt)
yline(0,'k--')
ylabel('$v(t)$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
grid on
legend({'$v_1(t)$','$v_2(t)$'},'Interpreter','latex','Location','best')
subplot(413)
plot(results.t_grid,I_opt)
ylabel('$I(t)$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
grid on
% results.t_grid_u = results.t_grid_u';
subplot(414)
results.u = [results.u,nan];
stairs(results.t_grid_u,results.u);
ylabel('$u(t)$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
grid on
% ylim([-1.1 1.1])
xlim([0 model.T])
