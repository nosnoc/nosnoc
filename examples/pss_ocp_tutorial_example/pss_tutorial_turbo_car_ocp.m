clear; clc; close all;
import casadi.*
import nosnoc.*
% 
problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
solver_options = nosnoc.solver.Options(); % Initialize all options related to the MPEC solver used for solving nosonc problems.

% Choosing the Runge - Kutta Method and number of stages
problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
problem_options.n_s = 2; % Number of stage points in the RK method (determines accuracy)

% Time-settings  - Solve a time optimal control problem
problem_options.time_optimal_problem = 1; 
problem_options.N_stages = 10; % Number of control intervals over the time horizon.
problem_options.N_finite_elements = 3; % Number of finite element (integration steps) on every control interval (optionally a vector might be passed).
problem_options.T = 1;    % Time horizon (for a time-optimal problem, the actual horizon length is obtained by solving the problem).

% settings.nlpsol = 'snopt';  % Note: requires installing.

model = nosnoc.model.Pss(); % Initialize a nosnoc model, (depending on the underlying nonsmooth model, several options exist)

% define differential states and populate the model.
q = SX.sym('q'); 
v = SX.sym('v'); % CasADi symbolic variables for states
model.x = [q;v]; % populate model state vectors
model.x0 = [0;0]; % initial value
v_max = 20;
model.lbx = [-inf;-v_max]; % lower bounds on states
model.ubx = [inf;v_max]; % upper bounds on states
% define control vectors
u = SX.sym('u');  % CasADi symbolic variables for controls
model.u = u;
u_max = 5;
model.lbu = -u_max ; 
model.ubu = u_max ;
% Dynamics of the piecewise smooth systems
f_1 = [v;u]; % mode 1 - nominal
f_2 = [v;3*u]; % mode 2 - turbo
% Define the regions of the PSS
v_threshold = 10;
model.c = v-v_threshold; % single switching functions (implies two regions)
model.S = [-1;1]; % Region R_1 is defined by c<0 (hence the -1), region R_2 by c>0 (hence the +1) in the sign matrix S.
model.F = [f_1 f_2]; % The columns of this matrix store the vector fields of every region.

model.g_terminal = [q-200;v-0]; % Add terminal constraint

% Create a nosnoc solver. 
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
% Solve the problem by internally calling the nosnoc MPEC solver (see generic_mpecs example for its standalone use)
ocp_solver.solve();
%% Extract reults - use ocp_solver methods to extact
t_grid = ocp_solver.get_time_grid(); % get time grid for differential states
t_grid_u = ocp_solver.get_control_grid(); % get time grid for control variables
x_opt = ocp_solver.get("x");  % get optimal solution for differential states
u_opt = ocp_solver.get("u"); % get optimal solutions for control variables
q_opt = x_opt(1,:);
v_opt = x_opt(2,:);

if problem_options.time_optimal_problem
    % Print value of optimal time.
    T = ocp_solver.get("T_final");
    fprintf('Final time: %2.4f s.\n',T)
else
    % In case the problem was not a time-optimal control problem, simply
    % print the objective.
    problem_options.T;
    fprintf('Objective value time: %2.4f s.\n',results.f)
end
%% Plot results
% enable latex formatting in plots
latexify_plot()
% Plot optimal states and controls
figure
subplot(131)
plot(t_grid,q_opt,'LineWidth',1.5);
xlabel('$t$','Interpreter','latex');
ylabel('$q(t)$','Interpreter','latex');
xlim([0 T])
grid on
subplot(132);
plot(t_grid,v_opt,'LineWidth',1.5);
yline(v_max,'r--');
yline(v_threshold,'k--');
ylim(1.05*[0 v_max]);
xlim([0 T])
xlabel('$t$','Interpreter','latex')
ylabel('$v(t)$','Interpreter','latex')
grid on
subplot(133);
stairs(t_grid_u,[u_opt,nan],'LineWidth',1.5);
xlabel('$t$','Interpreter','latex')
ylabel('$u(t)$','Interpreter','latex')
yline(-u_max,'r--');
yline(u_max,'r--');
ylim(1.2*[-u_max u_max]);
grid on
xlim([0 T])
