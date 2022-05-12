clear all
import casadi.*
[settings] = default_settings_fesd();  
% Choosing the Runge - Kutta Method and number of stages
settings.irk_scheme = 'Radau-IIA';
% settings.cross_comp_mode = 1;
settings.n_s = 2;
% Time-settings  - Solve an time optimal control problem
settings.time_optimal_problem = 1;
% Model - define all problem functions and
% Discretization parameters
model.N_stages = 10; % number of control intervals
model.N_finite_elements = 6; % number of finite element on every control intevral (optionally a vector might be passed)
model.T = 15;    % Time horizon
% Symbolic variables and bounds
q = SX.sym('q'); v = SX.sym('v'); 
model.x = [q;v]; % add all important data to the struct model,
model.x0 = [0;0]; % inital value
% bounds on states
model.lbx = [-inf;-20];
model.ubx = [inf;20];
% control
u = SX.sym('u'); model.u = u;
model.lbu = -5; model.ubu = 5;
% Dyanmics and the regions
f_1 = [v;u]; % mode 1 - nominal
f_2 = [v;3*u]; % mode 2 - turbo
model.S = [-1;1];
model.F = [f_1 f_2];
model.c = v-10;
% Add terminal constraint
model.g_terminal = [q-200;v-0];
% Solve OCP
[results,stats,model,settings] = nosnoc_solver(model,settings);

%% Plot
v_max = 20;
u_max = 5;
v_trash_hold = 10;
plot_results_nosnoc_tutorial
