clear; clc; close all;
import casadi.*
import nosnoc.*
%%
sliding_mode = 0;

problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
solver_options = mpecopt.Options(); % Initialize all options related to the MPEC solver used for solving nosonc problems.


problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
% problem_options.rk_scheme = RKSchemes.LOBATTO_IIIC;
problem_options.n_s = 2; % Number of stage points in the RK method (determines accuracy)

% Time-settings  - Solve a time optimal control problem
problem_options.time_optimal_problem = 0; 
problem_options.N_stages = 1; % Number of control intervals over the time horizon.
problem_options.N_finite_elements = 10; % Number of finite element (integration steps) on every control interval (optionally a vector might be passed).
problem_options.T = 2;    % Time horizon (for a time-optimal problem, the actual horizon length is obtained by solving the problem).
problem_options.relax_terminal_numerical_time = "ELL_1";
problem_options.rho_h = 1e-3;
problem_options.step_equilibration = "heuristic_diff";
% problem_options.rho_terminal_numerical_time = 1e4;
% problem_options.step_equilibration = "direct_homotopy_lift";
problem_options.cross_comp_mode = "FE_FE";
problem_options.dcs_mode = 'Heaviside';
% problem_options.cross_comp_mode = "STAGE_STAGE";

model = nosnoc.model.Pss(); % Initialize a nosnoc model, (depending on the underlying nonsmooth model, several options exist)
% define differential states and populate the model.
x = SX.sym('x'); 
model.x = x;
model.x0 = -1;
model.lbx = -inf; % lower bounds on states
model.ubx = inf;
% define control vectors
% Dynamics of the piecewise smooth systems
if sliding_mode 
    f_1 = 3;
    f_2 = -1;

    % active_set_guess = nosnoc.activeset.Heaviside({[1], [1 2]},'times', [problem_options.T/2 ,problem_options.T]);
    active_set_guess = nosnoc.activeset.Pss({[1], [1 2]},'times', [problem_options.T/2 ,problem_options.T]);

else
    f_1 = 3;
    f_2 = 1;
    % active_set_guess = nosnoc.activeset.Heaviside({[1], [2]},'times', [problem_options.T/2 ,problem_options.T]);
    active_set_guess = nosnoc.activeset.Pss({[1], [1 2]},'times', [problem_options.T/2 ,problem_options.T]);
    % active_set_guess = nosnoc.activeset.Heaviside({[1], [1 2]},'times', [problem_options.T/2 ,problem_options.T]);
end
% Define the regions of the PSS
model.c = x;
model.S = [-1;1]; % Region R_1 is defined by c<0 (hence the -1), region R_2 by c>0 (hence the +1) in the sign matrix S.
model.F = [f_1 f_2]; % The columns of this matrix store the vector fields of every region.


% Create a nosnoc solver.
solver_options.initialization_strategy = 'TakeProvidedActiveSet';
solver_options.rescale_large_objective_gradients = true;
solver_options.rho_TR_phase_ii_init = 1e-3;
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
% Set active set
ocp_solver.set_initial_active_set(active_set_guess);
ocp_solver.solve();


% Create a nosnoc solver using homotopy solver%ocp_solver = nosnoc.ocp.Solver(model, problem_options, homotopy_options);
% ocp_solver = nosnoc.ocp.Solver(model, problem_options, mpecopt_options);
% ocp_solver.set_initial_active_set(active_set_guess);
% ocp_solver.solve();
%% Extract reults - use ocp_solver methods to extact
t_grid = ocp_solver.get_time_grid(); % get time grid for differential states
x_opt = ocp_solver.get("x");  % get optimal solution for differential states

t_grid_full = ocp_solver.get_time_grid_full();
theta_opt = ocp_solver.get_full("theta");
lambda_opt = ocp_solver.get_full("lambda");

x_opt = x_opt(1,:);
figure
subplot(311)
plot(t_grid,x_opt);
xline(t_grid,'k')
subplot(312)
plot(t_grid_full,theta_opt,'.-')
hold on;
xline(t_grid,'k')
subplot(313)
plot(t_grid_full,lambda_opt,'.-')
hold on;
xline(t_grid,'k')
