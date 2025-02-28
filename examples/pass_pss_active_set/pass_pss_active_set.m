clear; clc; close all;
import casadi.*
import nosnoc.*
%% load options
problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
solver_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosonc problems.

%% set options
T = 25; % time horizon
N_stages = 10; % control stages
% Choosing the Runge - Kutta Method and number of stages
problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
problem_options.n_s = 1; % Number of stage points in the RK method (determines accuracy)

% Time-settings
problem_options.N_stages = N_stages; % Number of control intervals over the time horizon.
problem_options.N_finite_elements = 2; % Number of finite element (integration steps) on every control interval (optionally a vector might be passed).
problem_options.gamma_h = 1;
problem_options.T = T;    % Time horizon (for a time-optimal problem, the actual horizon length is obtained by solving the problem).

problem_options.relax_terminal_constraint = ConstraintRelaxationMode.ELL_1;
problem_options.cross_comp_mode = 'FE_FE';
problem_options.time_optimal_problem = 0;


% define differential states and populate the model.
q = SX.sym('q'); 
v = SX.sym('v'); 
u = SX.sym('u');  
v_max = 20;
u_max = 5;
% Dynamics of the piecewise smooth systems
f_1 = [v;u]; % mode 1 - nominal
f_2 = [v;3*u]; % mode 2 - turbo
% Define the regions of the PSS
v_threshold = 10;

model = nosnoc.model.Pss(); % Initialize a nosnoc model, (depending on the underlying nonsmooth model, several options exist)
model.x = [q;v]; % populate model state vectors
model.x0 = [0;0]; % initial value
model.lbx = [-inf;-v_max]; % lower bounds on states
model.ubx = [inf;v_max]; % upper bounds on states
model.u = u;
model.lbu = -u_max; 
model.ubu = u_max;
model.c = v-v_threshold; % single switching functions (implies two regions)
model.S = [-1;1]; % Region R_1 is defined by c<0 (hence the -1), region R_2 by c>0 (hence the +1) in the sign matrix S.
model.F = [f_1 f_2]; % The columns of this matrix store the vector fields of every region.
model.g_terminal = [q-200;v-0]; % Add terminal constraint
model.f_q = u^2;
%model.f_q_T = [q-200;v-0]'*diag([1000,1000])*[q-200;v-0]; % Add terminal constraint

% Define active set guess:
active_set_guess = nosnoc.activeset.Pss({[1],[2],[1]},'times', [T/4,3/4*T,T]);
active_set_guess = nosnoc.activeset.Pss({[1]},'times', [T]);
%active_set_guess = nosnoc.activeset.Pss({[1]},'times', [T-5]); % THROWS ERROR!
% Alternatively you can directly pass the indices of the switches as (control stage, finite element) pairs.
% This can be nice if you don't know anything about the specific times but know a general switching sequence.
% active_set_guess = nosnoc.activeset.Pss({[1],[2],[1]},'stages',{[3,2],[6,1],[10,3]});

% Create a nosnoc solver.
solver_options.mpecopt.initialization_strategy = 'TakeProvidedActiveSet';
solver_options.mpecopt.rescale_large_objective_gradients = true;
solver_options.mpecopt.rho_TR_phase_ii_init = 10;
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
% Set active set
ocp_solver.set_initial_active_set(active_set_guess);
ocp_solver.solve('mpecopt');

%% Extract reults - use ocp_solver methods to extact
t_grid = ocp_solver.get_time_grid(); % get time grid for differential states
t_grid_u = ocp_solver.get_control_grid(); % get time grid for control variables
x_opt = ocp_solver.get("x");  % get optimal solution for differential states
u_opt = ocp_solver.get("u"); % get optimal solutions for control variables
q_opt = x_opt(1,:);
v_opt = x_opt(2,:);

if problem_options.time_optimal_problem
    T = ocp_solver.get("T_final"); % Print value of optimal time.
    fprintf('Final time: %2.4f s.\n',T)
else
    % In case the problem was not a time-optimal control problem, simply print the objective.
    T = problem_options.T;
    fprintf('Objective value: %2.4f\n',ocp_solver.get_objective())
end
%% Plot results
% enable latex formatting in plots
latexify_plot()
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
