function [model] = car_hystheresis_model(model_in)

import casadi.*

unfold_struct(model_in,'caller');

%% Discretization parameters
N_stages = 30;
N_finite_elements = 3;


active_control = 1; % for ocps
% active_control = 0; % for sim

smooth_model = 0;
time_freezing_in_model = time_freezing;
%% goal
q_goal = 200;
v_goal = 0;
%% Bounds
v_max = 20;
umax = 5;

%% Time freezing and hysteresis settings
% normal ocp
if time_optimal_problem 
    N_stages = 35;
    N_finite_elements = 3;


    N_stages = 20;
    N_finite_elements = 3;
    % time optimal OcP
    fuel_cost_on = 0;
    control_cost_on = 0;
    v_trash_hold = 10; % switching criterion
    T = 3;   
else
    N_stages = 20;
    N_finite_elements = 3;
    fuel_cost_on = 0;
    control_cost_on = 1;
    v_trash_hold = 10; % switching criterion
    T = 15; % (here determined latter depeding on omega)
end


%% Model Parameters


% inital value
q0 = 0;
v0 = 0;
L0 =  0; % cost integral
a0 = 0;
t0 = 0;


% fual costs of turbo and nominal
% case 1
Pn = 1;
Pt = 1;
% case 2
Pn = 10;
Pt = 1;
%% Inital Value
x0 = [q0;v0;L0;a0;t0];

%% Time horizon
h = T/N_stages;

%% Model parameters for time freezing
% paramter of auxliary dynamcis
a_push = 5;
%% Define model dimensions, equations, constraint functions, regions an so on.
n_simplex = 1;% number of Carteisna products in the model ("independet switches"), we call this layer
% number of modes in every simplex
m_1 = 2;
m_vec = [m_1];
%% Variable defintion
q = MX.sym('q');
v = MX.sym('v');
L = MX.sym('L');
a = MX.sym('a');
t = MX.sym('t');

x = [q;v;L;a;t];
n_x = length(x);


if active_control
    lbx = -[inf;0;inf;inf;inf];
    ubx = [inf;v_max;inf;inf;inf];
else
    lbx = -inf*ones(n_x,1);
    ubx = inf*ones(n_x,1);
end
% every constraint funcion corresponds to a simplex (note that the c_i might be vector valued)


c_1 = v-v_trash_hold;


% sign matrix for the modes
S1 = [-1;...
       1];

% discrimnant functions
h_1 = -S1*[c_1];

c = [c_1];
h_indictaros = [h_1];

%% control
u = MX.sym('u');
n_u = 1;
u0 = [0];



if active_control
    lbu = -umax*ones(n_u,1);
    ubu = umax*ones(n_u,1);
else
    lbu = umax*ones(n_u,1);
    ubu = umax*ones(n_u,1);
end

%% modes of the ODEs layers (for all  i = 1,...,n_simplex);


if smooth_model
    f_nominal_1 = [v;3*u;Pn;0;1];
    f_nominal_2 = [v;3*u;Pn;0;1];
else
    f_nominal_1 = [v;u;Pn;1;1];
    f_nominal_2 = [v;3*u;Pt;0;1*(1-time_freezing_in_model)];
    f_nominal_2 = [v;3*u;Pt;0;1];
end

f_11 = f_nominal_1;
f_12 = f_nominal_2;

% in matrix form
f_1 = [f_11 f_12];

%% objective

f_q = fuel_cost_on*L+control_cost_on*u^2;
% Terminal Cost
f_q_T = 0;

%%  general nonlinear constinrst
general_nonlinear_constraint  = 0;
g_ineq = u^2;
g_ineq_lb = [-inf];
g_ineq_ub = [inf];
g_ineq_fun  = Function('g_ineq_fun',{x,u},{g_ineq});

%% terminal constraint
terminal_constraint = 1;

terminal_constraint =terminal_constraint*active_control;
g_terminal = [q-q_goal;v-v_goal];
g_terminal_lb = [0;0];
g_terminal_ub = [0;0];
g_terminal_fun = Function('g_terminal_fun',{x},{g_terminal});

%% Generic part
% (make of local workspace a struct and pass to output
names = who;
for ii = 1:length(names)
    eval([ 'model.' names{ii} '=' names{ii} ';'])
end
end

