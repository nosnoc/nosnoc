function [model] = car_hystheresis_model_voronoi(varargin)
import casadi.*
if nargin == 1
    unfold_struct(varargin{1},'caller')
else
    active_control = 1; % OCP or simulation problem
    linear_auxiliary_dynamics= 0; % constant or linear auxiliary dynamics
    smooth_model = 0;
    use_hystereis_model = 1;
    time_optimal_problem = 0;
    fuel_cost_on = 1;
    fuel_cost_same = 0;
end
fuel_cost_same = 0;
%% Terminal constraint and bounds
q_goal = 150;
v_goal = 0;
v_max = 25;
u_max = 5;
%% Discretization parameters
N_stages = 10;
N_finite_elements = 3;

T = 5;
%% Model Parameters
% inital value
q0 = 0;
v0 = 0;
L0 =  0; % cost integral
a0 = 0;
t0 = 0;
%% Hystheresis parameters
a1 = 0;
a2 = 1;
v1 = 10;
v2 = 15;
u0 = 0;
% fual costs of turbo and nominal
if fuel_cost_same
    Pn = 1;
    Pt = 1;
else
    Pn = 1;
    Pt = 10;
%     Pt = 100;
end
%% Inital Value
x0 = [q0;v0;L0;a0;t0];
%% Time horizon
h = T/N_stages;
%% Model parameters for time freezing

%% Define model dimensions, equations, constraint functions, regions an so on.
n_simplex = 1;% number of Carteisna products in the model ("independet switches"), we call this layer
if use_hystereis_model
    m_1 = 4;
else
    m_1 = 2;
end
m_vec = [m_1];
%% Variable defintion
q = MX.sym('q');
v = MX.sym('v');
L = MX.sym('L');
w = MX.sym('w');
t = MX.sym('t');


x = [q;v;L;w;t];

if active_control
    lbx = -[inf;0;inf;inf;inf];
    ubx = [inf;v_max;inf;inf;inf];
else
    lbx = -inf*ones(n_x,1);
    ubx = inf*ones(n_x,1);
end

k =  (a2-a1)/(v1-v2);
a_lin = a2-k*v1;
a_lin = a1-k*v2;

if ~use_hystereis_model
    c_1 = v-v2;
    S = [-1;1];  
end

%% PSS via Voronoi Cuts

z1 = [1/4;-1/4];
z2 = [1/4;1/4];
z3 = [3/4;3/4];
z4 = [3/4;5/4];
Z = [z1 z2 z3 z4];

psi = (v-v1)/(v2-v1);
z = [psi;w];

g_11 = norm([psi;w]-z1)^2;
g_12 = norm([psi;w]-z2)^2;
g_13 = norm([psi;w]-z3)^2;
g_14 = norm([psi;w]-z4)^2;
g_ind = [g_11;g_12;g_13;g_14];
g_ind_all = [g_ind];
c = g_ind_all;


%% control
u = MX.sym('u');
% u = SX.sym('u');


if active_control
    lbu = -u_max;
    ubu = u_max;
else
    lbu = u_max;
    ubu = u_max;
end
%% modes of the ODEs layers (for all  i = 1,...,n_simplex);
a_push_const = 0.5;
f_A = [v;u;Pn;0;1];
f_B = [v;3*u;Pt;0;1];
f_push_down = [0;0;0;-a_push_const;0];
f_push_up = [0;0;0;a_push_const;0];

if use_hystereis_model
    if smooth_model
        f_A = [v;u;Pn;0;1];
        f_11 = f_A;
        f_12 = f_A;
        f_13 = f_A;
        f_14 = f_A;
    else
            a_push = 1;
              f_push_down = [0;0;0;-a_push*(psi-1)^2/(1+(psi-1)^2);0];
              f_push_up = [0;0;0;a_push*(psi)^2/(1+(psi)^2);0];
              f_12 = f_push_down;
              f_13 = f_push_up;
              f_14 = 2*f_B-f_13;
              f_11 = 2*f_A-f_12;
    end
    % in matrix form
    f_1 = [f_11 f_12 f_13 f_14];
else
    if smooth_model
        f_11 = f_A;
        f_12 = f_A;
    else
        f_11 = f_A;
        f_12 = f_B;
    end
    f_1 = [f_11 f_12];
end
% matrix with all vector fields
F = f_1;
%% objective

% f_q = fuel_cost_on*L+(1-fuel_cost_on)*(1-time_optimal_problem)*u^2;
f_q = fuel_cost_on*L;
% Terminal Cost
% f_q_T = 0;

%%  general nonlinear constinrst
% no nonlinear constraints

%% terminal constraint
terminal_constraint = 1;
terminal_constraint =terminal_constraint*active_control;
g_terminal = [q-q_goal;v-v_goal];
% g_terminal_lb = [0;0];
% g_terminal_ub = [0;0];

%% Generic part
% (make of local workspace a struct and pass to output
names = who;
for ii = 1:length(names)
    eval([ 'model.' names{ii} '=' names{ii} ';'])
end
end

