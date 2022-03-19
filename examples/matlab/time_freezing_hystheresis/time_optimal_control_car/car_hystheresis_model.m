function [model] = car_hystheresis_model(varargin)
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
% N_stages = 30;
N_finite_elements = 3;
T = 12;

T = 5;
% T = 15;

%% Hystheresis parameters
a1 = 0;
a2 = 1;

v1 = 10;
v2 = 15;

% parameters of auxliary dynamics
a_push_const = 0.5;
a_push_linear = 0.25;
a_push_linear = 0.5/3;
a_push_linear = 0.5/5;

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
    Pt = 2;
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
    m_1 = 6;
else
    m_1 = 2;
end

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
k =  (a2-a1)/(v1-v2);
a_lin = a2-k*v1;
a_lin = a1-k*v2;




if use_hystereis_model
    c_1 = a-a1; % bottom
    c_2 = a-a2; % up
    c_3 = a-(k*v+a_lin); % diagonal
else
    c_1 = v-v2;
end


%%
if 0
    tt = -(2*v1):0.1:(2*v2);
    plot(tt,k*tt+a_lin)
    xline(v1)
    xline(v2)
    grid on
end

%%

% R1 - bottom with DAE forming for a = 0; c1 < 0 , c2 < 0 , c3 <0
% R2 - left, pushing down toward a = 0 ; c1 > 0 , c2 <0 ,c 3 < 0
% R3 - left, pushing down toward a = 0 ; c1 > 0 , c2 > 0 ,c 3 <0
% R4 - top with DAE forming for a = 1; c1 > 0 , c2 > 0 , c3 > 0
% R5 - right, pushing up towards a = 1; c1 > 0 ,c2 < 0, c3 > 0
% R6-  right, pushing up towards a = 1; c1 < 0 ,c2 > 0, c3 > 0

% sign matrix for the modes
if use_hystereis_model

    S1 = [-1 -1 -1;...
        1 -1 -1;...
        1 1 -1;...
        1 1 1;...
        1 -1 1;...
        -1 -1 1;...
        ];

    % discrimnant functions
    h_1 = -S1*[c_1;c_2;c_3];
    c = [c_1;c_2;c_2;c_3];
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
u0 = 0;

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
        f_15 = f_nominal_0;
        f_16 = f_nominal_0;
    else
        if linear_auxiliary_dynamics
            % for a = 0 and towards it, "left region"
            f_12 = [0;0;0;a_push_linear*(v-v2);0];
            f_13 = f_12;
            f_11 = 2*f_nominal_0-f_12;

            % for a = 1 and towards it, "rigt region"
            f_15 = [0;0;0;a_push_linear*(v-v1);0];
            f_16 = f_15 ;
            f_14 = 2*f_nominal_1-f_15 ;
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
    f_1 = [f_11 f_12 f_13 f_14 f_15 f_16];
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

