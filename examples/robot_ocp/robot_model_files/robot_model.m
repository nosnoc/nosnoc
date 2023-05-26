function [model] = robot_model()
import casadi.*

%% the dynamics
model.n_q = 4;
model.dims.n_dim_contact = 2;
g = -9.81*1;
% differential state
qx = MX.sym('qx',1);                                                       % position
qz = MX.sym('qz',1);                                                       % position
phi_hip = MX.sym('phi_hip',1);                                             % position
phi_knee = MX.sym('phi_knee',1);                                           % position
q = [qx;qz;phi_hip;phi_knee];
vx = MX.sym('vz',1);                                                       % velocity
vz = MX.sym('vz',1);                                                       % velocity
omega_hip = MX.sym('omega_hip',1);                                         % velocity
omega_knee = MX.sym('omega_knee',1);                                       % velocity
v = [vx;vz;omega_hip;omega_knee];
omega = v;
x = [q;v];
model.x = [q;v];
%% inital values
q0 = [0; 0.6 ;-pi/4;pi/2];
q0 = [0; 0.4 ;0;0.0];
v0 = [0;0;0.0;0.0];
model.x0 = [q0;v0];
model.u0 = [0;0];
%% friction cone parameters
model.e = 0.0;
model.mu_f = 0.7;
%%
mHip = 3.975; % mass of hip
mThigh = 1.782; %
mShank = 0.548;
%lengts
lBH = 0.043;  % distance between base and hip joint
lThigh = 0.2;
lShank = 0.2;
% center of masses distances
sBM = 0.02;  % distance between base and CoG of main body
sThigh = 0.016; % distance between hip joint and CoG of thigh *can be negative if z is positive?
sShank = 0.1;  % distance between knee joint and CoG of shank
rf = 0.028; % radius of foot
IyThugh = 0.001; % kgm2 inertia of thigh w.r.t. CoG about z-axis
IyShank = 0.0032;% inertia of shank w.r.t. CoG about z-axis
g = 9.81;
%% Control
u0 = [0;0];
tau_hip = MX.sym('tau_hip',1);                                             
tau_knee = MX.sym('tau_knee',1);                                             
model.u = [tau_hip;tau_knee];
control_force = [0;0;tau_hip;tau_knee];         
%% Dynamics and Kinematics
 % positions
% Hip
p_hipCOM = [q(1);q(2)];
% CoM
p_ThighCOM = p_hipCOM + sThigh*[-sin(q(3));-cos(q(3))]; 
p_ShankCOM = p_hipCOM + lThigh* [-sin(q(3));-cos(q(3))]+ sShank*[-sin(q(3)+q(4));-cos(q(3)+q(4))];
% Interesting Points
p_foot = p_hipCOM + lThigh* [-sin(q(3));-cos(q(3))]+ lShank*[-sin(q(3)+q(4));-cos(q(3)+q(4))]; 
p_knee = p_hipCOM + lThigh* [-sin(q(3));-cos(q(3))]; 
% velocity
v_hipCOM = p_hipCOM.jacobian(q)*v;
v_ThighCOM = jacobian(p_ThighCOM,q)*v;
v_ShankCOM = jacobian(p_ShankCOM,q)*v;
v_foot = jacobian(p_foot,q)*v;
v_knee = jacobian(p_knee,q)*v;
% casadi functins
p_foot_fun = Function('p_foot_fun',{x},{p_foot});
p_knee_fun = Function('p_knee_fun',{x},{p_knee});
v_foot_fun = Function('v_foot_fun',{x},{v_foot});
v_knee_fun = Function('v_knee_fun',{x},{v_knee});
v_hipCOM_fun  = Function('v_hipCOM_fun',{x},{v_hipCOM});
% intertia matrix
M_fun = Function.load('M_fun');          
M = M_fun([q;v]);
% bias forces
h_forces_fun= Function.load('h_forces_fun');
h_forces = h_forces_fun([q;v]);
% total forces unconstrained
model.f = h_forces+control_force;
%% constraint
model.f_c =  p_foot(2);
c_tan = p_foot(1);
model.tangent1 = c_tan.jacobian(q)';

end

