clear all
close all
import casadi.*

T = 20.0;
N_stages = 40;
N_finite_elements = 4;
%% Init nosnoc objects
model = nosnoc.model.PDSObjects; % Initialize model which is a container for the objects.
problem_options = nosnoc.Options; % Initialize all options related to the optimal control.
solver_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosonc problems.
%% DefinTe model
% Define controls and control bounds and add them to the model.
u = SX.sym('u', 2);
model.u = u;
model.lbu = [-1;-1];
model.ubu = [1;1];

ball1 = nosnoc.objects.Ball(0.25,2); % Define a ball in 2 dimensions with a radius of 0.25 units.

A1 = inv(diag([4, 0.1])); % Define the union of two ellipses with perpendicular minor axes of length 0.25
A2 = inv(diag([0.1, 4])); % and major axes of length 10.
ellipse2 = nosnoc.objects.Ellipse({A1,A2});

ball1.x0 = [-3;-3]; % Initial position of ball
ball1.f_rhs = u; % Ball is actuated directly with the controls.

ellipse2.x0 = [0;0;0]; % Initial position of the union of ellipses.
ellipse2.f_rhs = [0;0;0]; % Union of ellipses is unactuated.

model.addContact(ball1, ellipse2); % Add contact between ball and union of ellipses.

model.f_q = u'*diag([1e-1,1e-1])*u; % Add stage cost on controls
ball1_target = [-3;-3]; % define target state of the ball (its original position).
ellipse_target = [1;1;2*pi]; % Define target state of the union of ellipses, which is simply rotated one full rotation and shifted over by one unit in each direction.
x_target = vertcat(ball1_target, ellipse_target); % Total target state.
model.f_q_T = (model.x-x_target)'*diag([1e-1,1e-1,1e2,1e2,1e3])*(model.x-x_target); % Terminal cost penalizing deviaton from target state.

%% Define options
problem_options.T = T; % Time horizion
problem_options.N_stages = N_stages; % Number of control stages through the horizon.
problem_options.N_finite_elements = N_finite_elements; % Number of finite elements on each control stage.
problem_options.n_s = 2; % Number of steps in the Runge-Kutta Scheme.
problem_options.cross_comp_mode = CrossCompMode.FE_FE; % Fully dense cross complementarity mode
problem_options.gamma_h = 0.9; % Limit on finite element step variation which can aid in convergence at the possible cost of some optimality.


%% Define solver options.
solver_options.opts_casadi_nlp.ipopt.max_iter = 5000; % set max number of nlp iterations
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; % Use ma27 linear solver for performance (requires installadion)
solver_options.opts_casadi_nlp.ipopt.acceptable_iter = 5; % Only require three acceptable iter
solver_options.opts_casadi_nlp.ipopt.acceptable_tol = 1e-9; % Acceptable tolerance for nlp solver
solver_options.warm_start_duals = true; % Warm start nlp dual variables
solver_options.complementarity_tol = 1e-9; % Tolerance to solve complementarities to
%solver_options.polishing_step = 1;
solver_options.print_level = 3;

% Create a nosnoc solver. 
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
% Solve the problem by internally calling the nosnoc MPEC solver (see generic_mpecs example for its standalone use)
ocp_solver.solve();

x_res = ocp_solver.get('x'); % states
u_res = ocp_solver.get('u'); % controls 
pd_res = ocp_solver.get('p_d'); % contact point coordinates 
h_res = ocp_solver.get('h'); % step sizes

% Get matlab polygon objects from the nosnoc.objects types for plotting
pgon1 = ball1.to_polygon();
pgon2 = ellipse2.to_polygon();

% Define colors for plot
facecolor1 = [0 0.4470 0.7410];
linecolor1 = facecolor1*0.7;
facecolor2 = [0.8500 0.3250 0.0980];
linecolor2 = facecolor2*0.7;

% Create figure
fig = figure('Position', [10 10 1600 800]);

% Define indicies defining states.
indices = {1:2, 3:5};

% Plot example
plot_pds_sdf_example(h_res, x_res, pd_res, indices, [pgon1,pgon2], {facecolor1,facecolor2}, {linecolor1,linecolor2,}, fig, 'ellipse_spinning_friction');

