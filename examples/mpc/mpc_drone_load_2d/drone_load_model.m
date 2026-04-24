function [model, robot_data] = drone_load_model(params)
% DRONE_LOAD_MODEL  2D drone-load system with unilateral tether and ground contact.
%
% Description:
%   Planar aerial transport model for NOSNOC. A 2D drone is connected to a
%   passive load by a unilateral tether of length L. The tether can transmit
%   tension when taut, but no compression when slack. The load can also make
%   unilateral contact with the ground. No friction is included.
%
%   Contacts:
%     1) tether: active when ||p_d - p_l|| = L, transmits tension only
%     2) ground: active when the load touches the ground
%
%   States:
%       q = [x_d; z_d; phi; x_l; z_l]
%       v = [vx_d; vz_d; omega; vx_l; vz_l]
%       x = [q; v]
%
%     x_d, z_d : drone position [m]
%     phi      : drone pitch angle [rad]
%     x_l, z_l : load position [m]
%     vx_d, vz_d : drone velocity [m/s]
%     omega      : drone angular velocity [rad/s]
%     vx_l, vz_l : load velocity [m/s]
%
% Control:
%   If n_u = 2:
%       u = [T; tau]
%       T   : total thrust [N]
%       tau : pitch torque [Nm]
%
%   If n_u = 4:
%       u = [u1; u2; u3; u4]
%       individual rotor thrusts, internally mapped to total thrust and torque
%
%   The two/four inputs represent rotor thrusts. In the planar model they are
%   combined into total thrust and pitch torque acting on the drone.
%
%   Dynamics:
%       q_dot = v
%       M(q) v_dot = f_free(q,v,u) + J_n(q)' * lambda_n
%
%   The smooth dynamics contain drone thrust, torque, and gravity, and load
%   gravity. The contact forces consist of tether tension and ground normal
%   force.
%
%   The load is modeled as a point mass in the dynamics, while load_w and
%   load_h are used for plotting and for the ground contact offset.
if nargin < 1
    params = struct();
end

import casadi.*
model = nosnoc.model.Cls();

%% Parameters
g       = get_opt(params, 'g', 9.81);
m_d     = get_opt(params, 'm_d', 1.00);
I_d     = get_opt(params, 'I_d', 0.05);
m_l     = get_opt(params, 'm_l', 0.35);
L       = get_opt(params, 'L', 0.95);
arm     = get_opt(params, 'arm', 0.22);
u_max   = get_opt(params, 'u_max', 10.0);

c_d     = get_opt(params, 'c_d', 0.10);
c_phi   = get_opt(params, 'c_phi', 0.05);
c_l     = get_opt(params, 'c_l', 0.05);

drone_w = get_opt(params, 'drone_w', 0.40);
drone_h = get_opt(params, 'drone_h', 0.12);
load_w  = get_opt(params, 'load_w', 0.30);
load_h  = get_opt(params, 'load_h', 0.20);
ground_z = get_opt(params, 'ground_z', 0.0);

%% Dimensions
n_q = 5;
n_u = 2;
n_c = 2;
n_dim = 2;

%% Symbols
q = SX.sym('q', n_q);
v = SX.sym('v', n_q);
u = SX.sym('u', n_u);

%% Kinematics
p_d = [q(1); q(2)];
p_l = [q(4); q(5)];

%% Dynamics: M(q) * dv = f_v(q,v,u) + J_normal(q) * lambda_normal
% M: mass matrix
% f_v: smooth generalized forces = thrust + gravity + damping
% J_normal: contact Jacobian for tether and ground normal directions

M = SX.zeros(n_q, n_q);
M(1,1) = m_d;   % drone mass in x_d
M(2,2) = m_d;   % drone mass in z_d
M(3,3) = I_d;   % drone pitch inertia
M(4,4) = m_l;   % load mass in x_l
M(5,5) = m_l;   % load mass in z_l

if n_u == 4
    % four rotor thrusts, grouped into left/right thrust
    T_left  = u(1) + u(2);
    T_right = u(3) + u(4);
elseif n_u == 2
    % planar two-thruster model:
    % u(1) = left thrust, u(2) = right thrust
    T_left  = u(1);
    T_right = u(2);
else
    error('Only n_u = 2 or n_u = 4 supported.');
end

% total thrust and pitch torque
T_sum = T_left + T_right;
tau   = arm * (T_right - T_left);

% thrust in world frame
% phi = 0 => thrust points upward in +z direction
f_thrust = T_sum * [-sin(q(3)); cos(q(3))];

f_v = [ ...
    f_thrust(1) - c_d*v(1);         % drone x-force
    f_thrust(2) - m_d*g - c_d*v(2); % drone z-force
    tau - c_phi*v(3);               % drone pitch torque
   -c_l*v(4);                       % load x-force
   -m_l*g - c_l*v(5)];              % load z-force

%% Contact 1: unilateral tether (distance <= L)
r_dl = p_d - p_l;
dist_dl = sqrt(r_dl(1)^2 + r_dl(2)^2 + 1e-12);
e_dl = r_dl / dist_dl;

f_tether = L - dist_dl;
J_tether = [-e_dl(1);
            -e_dl(2);
             0;
             e_dl(1);
             e_dl(2)];

%% Contact 2: load-ground contact
f_ground = p_l(2) - ground_z - load_h/2;
J_ground = [0; 0; 0; 0; 1];

%% Initial state: load on the ground, tether taut, drone low
x_d0 = 0.0;
z_l0 = ground_z + load_h/2;
z_d0 = z_l0 + L;
q0 = [x_d0; z_d0; 0.0; 0.0; z_l0];
v0 = zeros(n_q,1);
x0 = [q0; v0];

%% Bounds
x_min = [-inf; -inf; -pi; -inf; -inf; 
         -40.0; -40.0; -100.0; -40.0; -40.0];

x_max = [ inf; inf; pi; inf; inf;
          40.0; 40.0; 100.0; 40.0; 40.0];

%% Pack NOSNOC model
model.q = q;
model.v = v;
model.u = u;
model.x = [q; v];
model.x0 = x0;

model.M = M;
model.f_v = f_v;
model.f_c = [f_tether; f_ground];
model.J_normal = [J_tether, J_ground];
model.J_tangent = SX.zeros(n_q, 0);
model.D_tangent = SX.zeros(n_q, 0);

model.mu = 0;
model.e = 0;
model.lbu = -u_max * ones(n_u,1);
model.ubu = u_max * ones(n_u,1);
model.lbx = x_min;
model.ubx = x_max;

%% Data for plotting / references
robot_data.g = g;
robot_data.m_d = m_d;
robot_data.I_d = I_d;
robot_data.m_l = m_l;
robot_data.L = L;
robot_data.arm = arm;
robot_data.u_max = u_max;
robot_data.ground_z = ground_z;
robot_data.drone_w = drone_w;
robot_data.drone_h = drone_h;
robot_data.load_w = load_w;
robot_data.load_h = load_h;
robot_data.n_q = n_q;
robot_data.n_u = n_u;
robot_data.n_c = n_c;
robot_data.n_dim = n_dim;
robot_data.x0 = x0;
robot_data.p_drone = p_d;
robot_data.p_load = p_l;
robot_data.hover_u = ((m_d + m_l) * g / 4) * ones(n_u,1);

end

function value = get_opt(s, name, default)
if isstruct(s) && isfield(s, name)
    value = s.(name);
else
    value = default;
end
end
