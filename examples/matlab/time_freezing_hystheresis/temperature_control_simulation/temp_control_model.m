function [model] = temp_control_model()

import casadi.*
%% Discretization parameters
N_stages = 2;
N_finite_elements = 1;
T = 0.1; % (here determined latter depeding on omega)
linear_push_dynamics = 1;
%% Model Parameters
active_control = 0;
% inital value
t0 = 0;
a = 0;
y0 = 12;

tau_cool_down = -0.3; % cool down time constant of lin dynamics
u_heat = 10; % heater power

% discrete values of the hysteresis function
a1 = 0;
a2 = 1;
% jump points in x in the hysteresis function
y1 = 18;
y2 = 20;

% slope and offset of middle diagonal line
k =  (a2-a1)/(y1-y2);
a_lin = a2-k*y1;
a_lin = a1-k*y2;

%% Inital Value

x0 = [y0;a;t0];

%% Time horizon
% N_stages = 75;
h = T/N_stages;

%% Model parameters for time freezing
a_push_const = 5;
a_push_linear = 1;


%% Define model dimensions, equations, constraint functions, regions an so on.
n_simplex = 1;% number of Carteisna products in the model ("independet switches"), we call this layer
% number of modes in every simplex
m_1 = 6;
m_vec = [m_1];
%% Variable defintion
y = MX.sym('y');
a = MX.sym('a');
t = MX.sym('t');

x = [y;a;t];
n_x = length(x);

lbx = -inf*ones(n_x,1);
ubx = inf*ones(n_x,1);
% every constraint funcion corresponds to a simplex (note that the c_i might be vector valued)




c_1 = a-a1; % bottom
c_2 = a-a2; % up
c_3 = a-k*y-a_lin; % diagonal

% R1 - bottom with DAE forming for a = 0; c1 < 0 , c2 < 0 , c3 <0
% R2 - left, pushing down toward a = 0 ; c1 > 0 , c2 <0 ,c 3 < 0
% R3 - left, pushing down toward a = 0 ; c1 > 0 , c2 > 0 ,c 3 <0
% R4 - top with DAE forming for a = 1; c1 > 0 , c2 > 0 , c3 > 0
% R5 - right, pushing up towards a = 1; c1 > 0 ,c2 < 0, c3 > 0
% R6-  right, pushing up towards a = 1; c1 < 0 ,c2 > 0, c3 > 0

% sign matrix for the modes
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
h_indictaros = [h_1];

%% control
u = MX.sym('u');
n_u = 1;
u0 = [0];

umax = 1e-3;

if active_control
    lbu = -umax*ones(n_u,1);
    ubu = umax*ones(n_u,1);
else
    lbu = 0*ones(n_u,1);
    ubu = 0*ones(n_u,1);
end

%% modes of the ODEs layers (for all  i = 1,...,n_simplex);
%
u_heat = 10;

f_nominal = [-0.2*y+(1-a)*u_heat;0;1];

f_nominal_A = [-0.2*y+u_heat;0;1];
f_nominal_B = [-0.2*y;0;1];


f_push_down = [0;-a_push_const;0];
f_push_up = [0;a_push_const;0];


a_gamma = 3;
f_push_down_gamma = [0;-a_gamma*(y-y2)^2/(1+(y-y2)^2);0];
f_push_up_gamma =  [0;a_gamma*(y-y1)^2/(1+(y-y1)^2);0];

if linear_push_dynamics
    % linear push dynamics
%     f_11 = 2*f_nominal+[0;-a_push_linear*(y-y2);0];
%     f_11 = 2*f_nominal_A+[0;-a_push_linear*(y-y2);0];
%     f_12 = [0;a_push_linear*(y-y2);0];
%     f_13 = [0;a_push_linear*(y-y2);0];
    f_11 = 2*f_nominal_A-f_push_down_gamma;
    f_12 = f_push_down_gamma;
    f_13 = f_push_down_gamma;
    % right toward a = 1
%     f_14 = 2*f_nominal+[0;-a_push_linear*(y-y1);0];
%     f_14 = 2*f_nominal_B+[0;-a_push_linear*(y-y1);0];
%     f_15 = [0;a_push_linear*(y-y1);0];
%     f_16 = [0;a_push_linear*(y-y1);0];
      f_14 = 2*f_nominal_B-f_push_up_gamma;
    f_15 = f_push_up_gamma;
    f_16 = f_push_up_gamma;
else
    % constant push dynamics
    % left toward a = 0
    f_11 = 2*f_nominal+f_push_up;
    f_12 = f_push_down;
    f_13 = f_push_down;
    % right toward a = 1
    f_14 = 2*f_nominal+f_push_down;
    f_15 = f_push_up;
    f_16 = f_push_up;
end
% in matrix form
f_1 = [f_11 f_12 f_13 f_14 f_15 f_16];

%% objective

f_q = active_control*(u^2)+y^2;

% Terminal Cost
f_q_T = 0;

%%  general nonlinear constinrst
general_nonlinear_constraint  = 0;
g_ineq = u^2;
g_ineq_lb = [-inf];
g_ineq_ub = [inf];
g_ineq_fun  = Function('g_ineq_fun',{x,u},{g_ineq});
%% Generic part
% (make of local workspace a struct and pass to output
names = who;
for ii = 1:length(names)
    eval([ 'model.' names{ii} '=' names{ii} ';'])
end
end

