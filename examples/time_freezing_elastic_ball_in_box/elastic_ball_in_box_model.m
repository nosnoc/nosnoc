function [model] = elastic_ball_in_box_model(omega)
import casadi.*
%% Discretization parameters
model = NosnocModel();
N_periods = 2;
alpha0 = pi/4; % inital angle
problem_options.T = N_periods*(2*pi/abs(omega));
%% Model Parameters
time_var_reference = 1;
qx_c = 0.0;
qy_c = 0.0;
v_target = 5;
R = 1;
friction_exists = 0;
u_max_R = 5e1; % amplitude of control force vector

% objective parameters
rho_q = 1;
rho_v = 1;
rho_u = 0;
%% Inital Value
qx0 = R*sin(alpha0);
qy0 = R*cos(alpha0);
vx0 = R*omega*cos(alpha0);
vy0 = -R*omega*sin(alpha0);
t0  = 0;
model.x0 = [qx0;qy0;vx0;vy0;t0];
%% Model parameters for time freezing
k_tf = 100;  % stiffnes 
gamma_tf = 1; % restitution coefficient
c_tf = 2*abs(log(gamma_tf))*((k_tf) /(pi^2+log(gamma_tf)^2) )^(1/2);
T_res = 2*pi/sqrt((4*k_tf-c_tf^2));
g = 9.81*0;
%% Variable defintion
qx = SX.sym('qx');
qy = SX.sym('qy');
vx = SX.sym('vx');
vy = SX.sym('vy');
t = SX.sym('t');
q = [qx;qy];
v = [vx;vy];
x = [q;v;t];
n_x = length(model.x);
n_q = 2;
model.lbx = -inf*ones(n_x,1);
model.ubx = inf*ones(n_x,1);

model.x = x;
%% control
ux = SX.sym('ux');
uy = SX.sym('uy');
model.u = [ux;uy];
u = model.u;
n_u = 2;
model.u0 = [0;0];
umax = inf;
model.lbu = -umax*ones(n_u,1);
model.ubu = umax*ones(n_u,1);
model.u = u;
%% Switching functions
% distance of constraints to (0,0) 
unit_size = 0.05*R*1;
b_bottom = -(R+1*unit_size);
a_right = (R+2*unit_size);
b_top = (R+3*unit_size);
a_left = -(R+4*unit_size);
% every constraint funcion corresponds to a simplex (note that the c_i might be vector valued)
c_1 = qy-b_bottom; % bottom
c_2 = -qx+a_right; % right
c_3 = -qy+b_top; % top
c_4 = qx-a_left; % left
% sign matrix for the modes
model.S = [1 1 1 1;...  % interior
     -1 1 1 1;...  % bottom
     -1 -1 1 1;...  % bottom right
      1 -1 1 1;...  % right
      1 -1 -1 1;...  % top right
      1 1 -1 1;...  % top
      1 1 -1 -1;...  % top left
      1 1 1 -1;...  % left
     -1 1 1 -1]; % bottom left
model.c = [c_1;c_2;c_3;c_4];

%% auxiliary dynamics
f_aux_right = [vx;0;-k_tf*(qx-a_right)-c_tf*vx;0;0];
f_aux_bottom = [0;vy;0;-k_tf*(qy-b_bottom)-c_tf*vy;0];
f_aux_left = [vx;0;-k_tf*(qx-a_left)-c_tf*vx;0;0];
f_aux_top = [0;vy;0;-k_tf*(qy-b_top)-c_tf*vy;0];
f_ode = [vx;vy;ux;uy-g;1];
f_zero = zeros(n_x,1);

f_11 = f_ode;
f_12 = f_aux_bottom;
f_13 = f_aux_bottom+f_aux_right;
f_14 = f_aux_right;
f_15 = f_aux_top+f_aux_right;
f_16 = f_aux_top;
f_17 = f_aux_top+f_aux_left;
f_18 = f_aux_left;
f_19 = f_aux_bottom+f_aux_left;
% in matrix form
model.F = [f_11 f_12 f_13 f_14 f_15 f_16 f_17 f_18 f_19];

%% objective
% if time_var_reference  
    qx_ref = R*sin(omega*t+alpha0);
    qy_ref = R*cos(omega*t+alpha0);    
    vx_ref = R*omega*cos(omega*t+alpha0);
    vy_ref = -R*omega*sin(omega*t+alpha0);
    q_ref = [qx_ref;qy_ref];    
    v_ref = [vx_ref;vy_ref];
    model.f_q = (rho_q*(q-q_ref)'*(q-q_ref)+rho_v*(v-v_ref)'*(v-v_ref)+rho_u*u'*u); 
% else
%     f_q = active_control*(rho_q*((qx-qx_c)^2+(qy-qy_c)^2-R^2)^2+rho_v*(v'*v-v_target^2)^2+rho_u*(u'*u));
% end
% Terminal Cost
model.f_q_T = 0;

%%  general nonlinear constinrst
model.g_path = u'*u;
model.g_path_lb = -inf;
model.g_path_ub = u_max_R^2;
end

