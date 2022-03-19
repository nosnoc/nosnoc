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

%% Terminal constraint and bounds
q_goal = 150;
v_goal = 0;
v_max = 25;
u_max = 5;
%% Discretization parameters
N_stages = 10;
N_finite_elements = 3;


% N_stages = 20;
% N_finite_elements = 3;


T = 5;

%% Hystheresis parameters
a1 = 0;
a2 = 1;

v1 = 10;
v2 = 15;
u0 = 0;


% parameters of auxliary dynamics
a_push_const = 0.5;
a_push_linear = 0.1;

%% Model Parameters
% inital value
q0 = 0;
v0 = 0;
L0 =  0; % cost integral
a0 = 0;
t0 = 0;

if 0
    % for simulation, "going down"
    u_max = -5;
    v0 = 25;
    a0 = 1;
end


% fual costs of turbo and nominal
if fuel_cost_same
    Pn = 1;
    Pt = 1;
else
    Pn = 1;
    Pt = 10;
end
%% Inital Value
x0 = [q0;v0;L0;a0;t0];

%% Time horizon
h = T/N_stages;

%% Model parameters for time freezing
% paramter of auxliary dynamcis

%% Define model dimensions, equations, constraint functions, regions an so on.
n_simplex = 1;% number of Carteisna products in the model ("independet switches"), we call this layer
% number of modes in every simplex
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
n_x = length(x);


if active_control
    lbx = -[inf;0;inf;inf;inf];
    ubx = [inf;v_max;inf;inf;inf];
else
    lbx = -inf*ones(n_x,1);
    ubx = inf*ones(n_x,1);
end
% every constraint funcion corresponds to a simplex (note that the c_i might be vector valued)
k =  (a2-a1)/(v1-v2);
a_lin = a2-k*v1;
a_lin = a1-k*v2;




if use_hystereis_model
%     c_1 = a-a1; % bottom
%     c_2 = a-a2; % up
%     c_3 = a-(k*v+a_lin); % diagonal
else
    c_1 = v-v2;
end

%% PSS via Voronoi Cuts

z1 = [1/4;-1/4];
z2 = [1/4;1/4];
z3 = [3/4;3/4];
z4 = [3/4;5/4];

Z = [z1 z2 z3 z4];

% sign matrix for the modes
if use_hystereis_model

psi = (v-v1)/(v2-v1);
z = [psi;w];
% discriminant functions via voronoi
h_1 = -2*z'*z1+norm(z1)^2;
h_2 = -2*z'*z2+norm(z2)^2;
h_3 = -2*z'*z3+norm(z3)^2;
h_4 = -2*z'*z4+norm(z4)^2;

% 
h_11 = norm([psi;w]-z1)^2;
h_12 = norm([psi;w]-z2)^2;
h_13 = norm([psi;w]-z3)^2;
h_14 = norm([psi;w]-z4)^2;

% h_1 = 1;
% h_2 = 1;
% h_3 = 15;
% h_4 = 25;

h_1 = [h_11;h_12;h_13;h_14];
h_indictaros = [h_1];
c = h_indictaros;
else
    S1 = [-1; ...
        1];
    % discrimnant functions
    h_1 = -S1*[c_1];
    c = [c_1];

end

h_indictaros = [h_1];

%% control
u = MX.sym('u');
n_u = 1;

if active_control
    lbu = -u_max*ones(n_u,1);
    ubu = u_max*ones(n_u,1);
else
    lbu = u_max*ones(n_u,1);
    ubu = u_max*ones(n_u,1);
end

%% modes of the ODEs layers (for all  i = 1,...,n_simplex);


f_nominal_0 = [v;u;Pn;0;1];
f_nominal_1 = [v;3*u;Pt;0;1];
f_push_down = [0;0;0;-a_push_const;0];
f_push_up = [0;0;0;a_push_const;0];


if use_hystereis_model
    if smooth_model
        f_nominal_0 = [v;u;Pn;0;1];
        f_11 = f_nominal_0;
        f_12 = f_nominal_0;
        f_13 = f_nominal_0;
        f_14 = f_nominal_0;
    else
        if linear_auxiliary_dynamics
%             % for a = 0 and towards it, "left region"
%             f_12 = [0;0;0;a_push_linear*(v-v2);0];
%             
% 
%             % for a = 1 and towards it, "rigt region"
%             f_13 = [0;0;0;a_push_linear*(v-v1);0];
              
              a_push = 1;
              f_push_down = [0;0;0;-a_push*(psi-1)^2/(1+(psi-1)^2);0];
              f_push_up = [0;0;0;a_push*(psi)^2/(1+(psi)^2);0];
              f_12 = f_push_down;
              f_13 = f_push_up;
              f_14 = 2*f_nominal_1-f_13;
              f_11 = 2*f_nominal_0-f_12;
        else
            % constant push dynamics
            % left toward a = 0
            f_11 = 2*f_nominal_0+f_push_up;
            f_12 = f_push_down;
            f_13 = f_push_down;
            % right toward a = 1
            f_14 = 2*f_nominal_1+f_push_down;
            f_15 = f_push_up;
            f_16 = f_push_up;
        end
    end


    % in matrix form
    f_1 = [f_11 f_12 f_13 f_14];
else
    if smooth_model
        f_11 = f_nominal_0;
        f_12 = f_nominal_0;
    else
        f_11 = f_nominal_0;
        f_12 = f_nominal_1;
    end
    f_1 = [f_11 f_12];

end
%% objective

f_q = fuel_cost_on*L+(1-fuel_cost_on)*(1-time_optimal_problem)*u^2;
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

