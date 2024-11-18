function [model] = car_hysteresis_model_voronoi(varargin)
import casadi.*
if nargin == 1
    unfold_struct(varargin{1},'caller')
else
    fuel_cost_on = 0;
    fuel_cost_same = 0;
end
model = nosnoc.model.Pss();
%% Terminal constraint and bounds
q_goal = 150;
v_goal = 0;
v_max = 25;
u_max = 5;
%% Model Parameters
% Hystheresis parameters
v1 = 10;
v2 = 15;
% fual costs of turbo and nominal
Pn = 1;
Pt = 1;

%% Variable defintion
% states
q = SX.sym('q');
v = SX.sym('v');
L = SX.sym('L');
w = SX.sym('w');
t = SX.sym('t');
x = [q;v;L;w;t];
% controls
u= SX.sym('u');

model.x = x;
model.u = u;

% Bounds on x and u
model.lbx = -[inf;0;inf;inf;inf];
model.ubx = [inf;v_max;inf;inf;inf];
model.lbu = -u_max;
model.ubu = u_max;
model.u0 = 4;
%% Inital Value
model.x0 = zeros(5,1);
% u0 = 10;
%% PSS via Voronoi Cuts
z1 = [1/4;-1/4];
z2 = [1/4;1/4];
z3 = [3/4;3/4];
z4 = [3/4;5/4];

psi = (v-v1)/(v2-v1);

g_1 = norm([psi;w]-z1)^2;
g_2 = norm([psi;w]-z2)^2;
g_3 = norm([psi;w]-z3)^2;
g_4 = norm([psi;w]-z4)^2;

model.g_indicator = [g_1;g_2;g_3;g_4];

% modes of the ODEs layers
f_A = [v;u;Pn;0;1];
f_B = [v;3*u;Pt;0;1];

a_push = 1;
gamma_p = a_push*(psi^2/(1+psi^2))
gamma_n = a_push*((psi-1)^2/(1+(psi-1)^2))

f_push_down = [0;0;0;-gamma_n;0];
f_push_up = [0;0;0;gamma_p;0];

f_2 = f_push_down;
f_3 = f_push_up;
f_4 = 2*f_B-f_3;
f_1 = 2*f_A-f_2;
% in matrix form
model.F = [f_1 f_2 f_3 f_4];
%% objective and terminal constraint
model.f_q = fuel_cost_on*L;
% terminal constraint
model.g_terminal = [q-q_goal;v-v_goal];
%model.a_n = 1000;
% g_terminal_lb = zeros(2,1);
% g_terminal_ub = zeros(2,1);
% f_q_T = 1e3*[q-q_goal;v-v_goal]'*[q-q_goal;v-v_goal];
end

