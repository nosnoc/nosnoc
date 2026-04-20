function [model, robot_data] = spinner_model(shape)
%--------------------------------------------------------------------------
% 2-Link planar spinner finger + rotating sphere contact
% Compatible with MATLAB NOSNOC (CasADi-based)
%--------------------------------------------------------------------------

if nargin < 1
    shape = 'circle';
end

import casadi.*
model = nosnoc.model.Cls();

%% Parameters
g  = 9.81;
mu = 0.8;           % friction coefficient

% Geometry
l1 = 1.0;           % upper link length
l2 = 1.0;           % lower link length
cx = 1.5;
cz = 1.0; % ball center (fixed in world)

switch lower(shape)
    case 'circle'
        R = 0.3;
        geom_type = 'circle';
        a = R; b = R;
    case 'ellipse'
        a = 0.3; b = 0.2;     % semi-axes
        R = a;
        geom_type = 'ellipse';
    otherwise
        error('Unknown shape type.');
end



% Mass/inertia (approx. from URDF)
m1 = 0.5;  I1 = 0.0839;
m2 = 0.5;  I2 = 0.0839;
mb = 0.5;  Ib = 0.074;

if strcmpi(shape,'circle')
    Ib = 0.5 * mb * a^2;          % classic disk
else
    Ib = 0.25 * mb * (a^2 + b^2); % ellipse
end


% COM offsets
c1 = l1/2;
c2 = l2/2;

%% Initial configuration
% fingertip target = leftmost point of shape (world x direction)
xt = cx - a;   % leftmost x
zt = cz;       % same z

% compute IK for 2-link arm
r2 = xt^2 + zt^2;  r = sqrt(r2);
r_min = abs(l1 - l2) + 1e-6;
r_max = (l1 + l2) - 1e-6;
r_clamp = min(max(r, r_min), r_max);

c2 = (r_clamp^2 - l1^2 - l2^2)/(2*l1*l2);
c2 = max(-1, min(1, c2));
q2_ik = acos(c2);                  % elbow-down
phi = atan2(zt, xt);
psi = atan2(l2*sin(q2_ik), l1 + l2*c2);
q1_ik = phi - psi;

q0 = [q1_ik; q2_ik; 0];
v0 = zeros(3,1);
x0 = [q0; v0];
model.x0 = x0;

% store base geometry for robot_data later
robot_data.a = a;  robot_data.b = b;
robot_data.cxy = [cx; cz];
robot_data.shape = geom_type;
robot_data.x0 = x0;


%% Dimensions
n_q = 3;     % [q1; q2; q3] -> shoulder, elbow, ball
n_u = 2;     % torques on shoulder and elbow
n_dim = 2;   % planar
n_c = 1;     % one contact point

%% Bounds

x_max = [ 2*pi; 2*pi; 10*pi; 8; 8; 8];

u_max = [20; 20];
x_max = [ 2*pi; 2*pi; 10*pi; 8; 8; 8];

u_max = [10; 10];
x_max = [ 10*pi; 10*pi; inf; inf; inf; inf];

x_min = -x_max;

%% CasADi symbols
q = SX.sym('q', n_q);
v = SX.sym('v', n_q);
u = SX.sym('u', n_u);

%% Kinematics -------------------------------------------------------------
p0 = [0; 0];
p1 = p0 + [l1*sin(q(1)); l1*cos(q(1))];
p2 = p1 + [l2*sin(q(1)+q(2)); l2*cos(q(1)+q(2))];   % fingertip
p_ball = [cx; cz];

% COM positions
p_c1 = p0 + [c1*sin(q(1)); c1*cos(q(1))];
p_c2 = p1 + [c2*sin(q(1)+q(2)); c2*cos(q(1)+q(2))];

% Jacobians of COMs
J_c1 = [ c1*cos(q(1)), 0, 0;
    -c1*sin(q(1)), 0, 0];
J_c2 = [ l1*cos(q(1))+c2*cos(q(1)+q(2)),  c2*cos(q(1)+q(2)), 0;
    -l1*sin(q(1))-c2*sin(q(1)+q(2)), -c2*sin(q(1)+q(2)), 0];

% Tip Jacobian
J_tip = [ l1*cos(q(1))+l2*cos(q(1)+q(2)),  l2*cos(q(1)+q(2)), 0;
    -l1*sin(q(1))-l2*sin(q(1)+q(2)), -l2*sin(q(1)+q(2)), 0];

%% Lagrangian dynamics ----------------------------------------------------
% Velocities
v_c1 = J_c1*v;
v_c2 = J_c2*v;

% Kinetic + potential energy
T = 0.5*m1*(v_c1.'*v_c1) + 0.5*m2*(v_c2.'*v_c2) ...
    + 0.5*I1*v(1)^2 + 0.5*I2*(v(1)+v(2))^2 + 0.5*Ib*v(3)^2;
V = m1*g*p_c1(2) + m2*g*p_c2(2);
L = T - V;

dLdq  = jacobian(L,q);
dLddq = jacobian(L,v);

M = jacobian(dLddq,v);                 % mass matrix
C = jacobian(dLddq,q)*v - dLdq.';      % Coriolis + gravity term

%% Inputs -----------------------------------------------------------------
B = [1 0 0;
    0 1 0];
% f_v = -C + B.'*u;

b_spinner = 0.1;  % damping coefficient (tune as desired)
C_fric = [0; 0; b_spinner*v(3)];
f_v = -C - C_fric + B.'*u;

%% Contact (tip ↔ sphere) -----------------------------------------------
% r = p2 - p_ball;                   % vector from ball center to fingertip
% r_norm = sqrt(r(1)^2 + r(2)^2);
% f_c = r_norm - R;                  % scalar gap
%
% % n = r / (r_norm + 1e-12);          % normal (outward from ball)
% % t = [-n(2); n(1)];                 % tangent
% %
% % % Contact Jacobians (n_q × n_c)
% % % J_normal  = [ J_tip.' * n ; SX(0) ];          % (3×1)
% % % J_tangent = [ J_tip.' * t ; SX(R) ];          % (3×1)
% % %
% % % J_normal  = [ J_tip.' * n];          % (3×1)
% % % J_tangent = [ J_tip.' * t];          % (3×1)
% %
% % J_normal  = J_tip.' * n;                      % (3×1)
% % J_tangent = J_tip.' * t;                      % (3×1)
% % J_tangent(3) = R;                             % add torque on spinner
%
%
% % unit normal/tangent at contact
% n = (p2 - p_ball) / (sqrt((p2(1)-p_ball(1))^2 + (p2(2)-p_ball(2))^2) + 1e-12);
% t = [-n(2); n(1)];                  % CCW +90° from n
%
% J_normal   = J_tip.'*n;             % (3x1)
% J_tangent  = J_tip.'*t;             % (3x1)
% J_tangent(3) = -R;                  % <-- flip sign here
%

if 0
    p_rel = p2 - [cx; cz];

    switch geom_type
        case 'circle'
            % implicit F = (x/a)^2 + (z/b)^2 - 1
            r_norm = sqrt(p_rel(1)^2 + p_rel(2)^2);
            f_c = r_norm - a;
            n = p_rel / (r_norm + 1e-12);

        case 'ellipse'
            F = (p_rel(1)/a)^2 + (p_rel(2)/b)^2 - 1;
            f_c = F;  % scalar gap >0: outside
            n = [2*p_rel(1)/(a^2); 2*p_rel(2)/(b^2)];
            n = n / (sqrt(n(1)^2 + n(2)^2) + 1e-12);
    end

    t = [-n(2); n(1)];          % tangent (+90°)
    J_normal  = J_tip.'*n;
    J_tangent = J_tip.'*t;

    % torque coupling on spinner
    r_eq = sqrt(a*b);            % effective radius
    J_tangent(3) = -r_eq;        % torque sign convention

else
    % --- Contact geometry: circle uses simple world formula; ellipse uses rotating local frame ---
    p_rel_world = p2 - [cx; cz];
    theta = q(3);

    if abs(a-b) < 1e-12
        % ===== CIRCLE (a=b=R): simple, rotation-invariant =====
        Rloc = a;                                   % radius
        r_norm = sqrt(p_rel_world(1)^2 + p_rel_world(2)^2);
        f_c = r_norm - Rloc;                        % gap
        n = p_rel_world / (r_norm + 1e-12);         % outward unit normal in world
    else
        % ===== ELLIPSE (a≠b): rotate world->local, compute F, map gradient back =====
        R_q = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        p_local = R_q.' * p_rel_world;              % world -> spinner-local

        % implicit function in local coords
        F = (p_local(1)/a)^2 + (p_local(2)/b)^2 - 1;
        f_c = F;

        gradF_local = [ 2*p_local(1)/(a^2) ;
            2*p_local(2)/(b^2) ];
        n = R_q * gradF_local;                      % back to world
        n = n / (sqrt(n(1)^2 + n(2)^2) + 1e-12);    % normalize
    end

    % tangent (+90° in world)
    t = [-n(2); n(1)];

    % contact Jacobians
    J_normal  = J_tip.' * n;                        % (3×1)
    J_tangent = J_tip.' * t;                        % (3×1)

    % spinner torque coupling from tangential force λ_t:
    %   τ3 = -r_eq * λ_t
    if abs(a-b) < 1e-12
        r_eq = a;                                   % circle: exact
    else
        r_eq = sqrt(a*b);                           % ellipse: effective radius
    end
    J_tangent(3) = -r_eq;
end


%% Initial configuration --------------------------------------------------
% xt = cx - R;  zt = cz;
% r2 = xt^2 + zt^2;  rlen = sqrt(r2);
% c2ang = (r2 - l1^2 - l2^2)/(2*l1*l2); c2ang = max(min(c2ang,1),-1);
% q2_ik = acos(c2ang);
% q1_ik = atan2(zt,xt) - atan2(l2*sin(q2_ik), l1 + l2*c2ang);
% q0 = [q1_ik; q2_ik; 0];
% v0 = zeros(n_q,1);
% x0 = [q0; v0];

%% Pack model -------------------------------------------------------------
model.q = q; model.v = v; model.u = u;
model.x = [q; v];
model.x0 = x0;

model.M = M;
model.f_v = f_v;
model.f_c = f_c;
model.J_normal  = J_normal;
model.J_tangent = J_tangent;

model.mu = mu;
model.e = 0;
model.lbu = -u_max; model.ubu = u_max;
model.lbx = x_min;  model.ubx = x_max;

% crucial: specify dimensions
% model.n_dim = n_dim;
% model.n_c   = n_c;

%% Robot data structure
robot_data.g = g; robot_data.mu = mu;
robot_data.l1 = l1; robot_data.l2 = l2;
robot_data.c1 = c1; robot_data.c2 = c2;
robot_data.p2 = p2;
robot_data.R = R;
robot_data.shape = geom_type;
robot_data.a = a;
robot_data.b = b;
robot_data.cxy = [cx; cz];
robot_data.m = [m1; m2; mb];
robot_data.I = [I1; I2; Ib];
robot_data.n_q = n_q;
robot_data.n_u = n_u;
robot_data.n_dim = n_dim;
robot_data.n_c = n_c;
robot_data.x0 = x0;

end
