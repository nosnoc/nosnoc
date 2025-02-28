clear; clc; close all;
import casadi.*
import nosnoc.*
% 
problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
solver_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosonc problems.

% Choosing the Runge - Kutta Method and number of stages
problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
problem_options.n_s = 2; % Number of stage points in the RK method (determines accuracy)

% Time-settings  - Solve a time optimal control problem
problem_options.time_optimal_problem = 1; 
problem_options.N_stages = 10; % Number of control intervals over the time horizon.
problem_options.N_finite_elements = 3; % Number of finite element (integration steps) on every control interval (optionally a vector might be passed).
problem_options.T = 1;    % Time horizon (for a time-optimal problem, the actual horizon length is obtained by solving the problem).

solver_options.print_level = 3;

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
% Generate the relevant c++ files as well as a CMakeLists.txt which expects a main.cpp provided by the user.
% The argument here generates the files in the current working directory, though this can be re-directed to any directory.
%
% Note: this relies on Matlab-Python interop. If that fails you will need to call the python script in `codegen_templates`
%       manually to generate the c++ files from the json information.
ocp_solver.generate_cpp_solver(fullfile(pwd));

% Now the focus switches to the c++ code, see main.cpp for details on how to use the generated solver.
% In order to build the c++ solver:
%
% 1. in matlab shell:
% > cpp_solver_tutorial
%
% 2. in system shell:
% > export CASADI_CMAKE_PATH=<path/to/casadi>/cmake 
% > mkdir build; cd build
% > cmake ..
% > make


