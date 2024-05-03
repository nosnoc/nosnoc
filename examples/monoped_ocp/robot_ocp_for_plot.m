import casadi.*
clear all;
close all;
%% robot scene description
constant_inertia_matrix = 0;
general_inequality_constraints = 0;
save_figure = 1;
filename = 'monoped_ocp';

%% auxiliary dynamics and friction
a_n = 200;
mu = 0.80;
%% obstacles
q_target = [3;0.4;0;0];

%% Default settings NOSNOC
settings = NosnocOptions();
model = NosnocModel();
%%
settings.print_level = 5;
settings.irk_scheme = IRKSchemes.RADAU_IIA;
settings.n_s = 2;
%% homotopy settings
settings.cross_comp_mode = 1;
settings.opts_casadi_nlp.ipopt.max_iter = 10000;
%settings.opts_casadi_nlp.ipopt.max_iter = 1000;
settings.N_homotopy = 10;
settings.sigma_0 = 1;
% settings.homotopy_update_rule = 'superlinear';
settings.homotopy_update_slope = 0.4;
settings.opts_casadi_nlp.ipopt.tol = 1e-5;
settings.opts_casadi_nlp.ipopt.acceptable_tol = 1e-5;
settings.opts_casadi_nlp.ipopt.acceptable_iter = 3;
settings.complementarity_tol = 1e-5;
settings.opts_casadi_nlp.ipopt.linear_solver = 'ma97';
% settings.nlpsol = 'snopt';
% settings.opts_casadi_nlp.snopt.Major_iterations_limit = 100000;
% settings.opts_casadi_nlp.snopt.Minor_iterations_limit = 10000;
settings.initial_alpha = 0.5;
%settings.opts_casadi_nlp.ipopt.ma57_pre_alloc = 2.0;
%settings.opts_casadi_nlp.ipopt.ma57_automatic_scaling = 'yes';

%% time-freezing
settings.s_sot_max = 100;
settings.s_sot_min = 0.99;
settings.rho_sot = 0.00;
settings.time_freezing = 1;
settings.pss_lift_step_functions = 0;
settings.stagewise_clock_constraint = 1;

%% Discretization
problem_options.T = 3.0;
settings.N_stages = 60;
settings.N_finite_elements = 3;

%% friction cone parameters
model.e = 0;
model.mu_f = mu;
model.a_n = a_n;
%% bounds
lb_head_z = 0.2;
ub_head_z = 0.55;
lb_head_x = -0.05;
ub_head_x = inf;
%
lb_knee_x = -0.05;
ub_knee_x = inf;
lb_knee_z = 0.05;
ub_knee_z = inf;
%
lb_foot_x = -0.05;
ub_foot_x = inf;
lb_foot_z = -0.005;
ub_foot_z = 0.2;

psi_hip_ub = 3*pi/8*1.05;
psi_hip_lb = -3*pi/8*1.05;
psi_knee_ub = pi/2*1.05;
psi_knee_lb = -pi/2*1.05;

%% robot model parameters
mHip = 3.975; % mass of hip
mThigh = 1.782; %
mShank = 0.548;
%lengts
lBH = 0.043;  % distance between base and hip joint
lThigh = 0.2;
lShank = 0.2;
lHead = 0.05;
% center of masses distances
sBM = 0.02;  % distance between base and CoG of main body
sThigh = 0.016; % distance between hip joint and CoG of thigh *can be negative if z is positive?
sShank = 0.1;  % distance between knee joint and CoG of shank
rf = 0.028; % radius of foot
IyThugh = 0.001; % kgm2 inertia of thigh w.r.t. CoG about z-axis
IyShank = 0.0032;% inertia of shank w.r.t. CoG about z-axis
g = 9.81;
%% differential state
qx = SX.sym('qx',1);
qz = SX.sym('qz',1);
phi_hip = SX.sym('phi_hip',1);
phi_knee = SX.sym('phi_knee',1);
vx = SX.sym('vz',1);
vz = SX.sym('vz',1);
omega_hip = SX.sym('omega_hip',1);
omega_knee = SX.sym('omega_knee',1);
% controls
u_hip = SX.sym('u_hip',1);
u_knee = SX.sym('u_knee',1);

u = [u_hip;u_knee];
q = [qx;qz;phi_hip;phi_knee];
v = [vx;vz;omega_hip;omega_knee];
x = [q;v];
model.x = [q;v];
model.q = q;
model.v = v;
model.u = u;
%% inital values
q0 = [0;0.4;0;0];
v0 = [0;0;0;0];
model.x0 = [q0;v0];
%% Dynamics and Kinematics
robot_model_kinematics
% total forces unconstrained
if constant_inertia_matrix
    q_lin = [0;0.4;pi/2;-pi/4];
    M = full(M_fun([q_lin;v0]));
end
model.f_v = (h_forces+[0;0;u]);
model.invM = inv(M);
model.M = M;
%% normal and tangents
f_c =  p_foot(2);
c_tan = p_foot(1);
J_normal = f_c.jacobian(q)';
J_tangent = c_tan.jacobian(q)';
use_unit_vectors = 0;
if use_unit_vectors
    J_tangent = J_tangent/norm(J_tangent);
    J_normal = J_normal/norm(J_normal);
end
model.f_c = f_c;
model.J_tangent = J_tangent;
model.J_normal= J_normal;
model.dims.n_dim_contact = 2;
%% OCP
% Objective and constraints
% box constraints
u_max = 100;
model.lbu = -u_max*ones(2,1);
model.ubu = u_max*ones(2,1);
% Sanity constraints
model.lbx = [-0.5;0;-pi;-pi;-100*ones(4,1)];
model.ubx = [q_target(1)+0.5; 10;pi;pi;100*ones(4,1)];
%% path constraints
% lower bound on knee
p_knee_x = p_knee(1);
p_knee_z = p_knee(2);

p_foot_x = p_foot(1);
p_foot_z = p_foot(2);

g_path = [];
g_path_lb = [];
g_path_ub = [];
% constraint on knee x
if (lb_knee_x ~= -inf) || (ub_knee_x ~= inf)
    g_path = [g_path;p_knee_x];
    g_path_lb = [g_path_lb;lb_knee_x];
    g_path_ub = [g_path_ub;ub_knee_x];
end
% constraint on knee z
if (lb_knee_z ~= -inf) || (ub_knee_z ~= inf)
    g_path = [g_path;p_knee_z];
    g_path_lb = [g_path_lb;lb_knee_z];
    g_path_ub = [g_path_ub;ub_knee_z];
end

% constraint on foot x
if (lb_foot_x ~= -inf) || (ub_foot_x ~= inf)
    g_path = [g_path;p_foot_x];
    g_path_lb = [g_path_lb;lb_foot_x];
    g_path_ub = [g_path_ub;ub_foot_x];
end

% constraint on foot z
if (lb_foot_z ~= -inf) || (ub_foot_z ~= inf)
    g_path = [g_path;p_foot_z];
    g_path_lb = [g_path_lb;lb_foot_z];
    g_path_ub = [g_path_ub;ub_foot_z];
end

if general_inequality_constraints
    if ~isempty(g_path)
        model.g_path = g_path;
        model.g_path_ub = g_path_ub;
        model.g_path_lb = g_path_lb;
    end
end

% least squares weight

Q = diag([1, 1, 10, 1, 1e-6, 1e-6, 1e-6, 1e-6]);
Q_terminal = diag([1e3, 1e3, 1e3, 1e3, 10, 10, 10, 10]);


u_ref = [0;0];

R = 1e-1*eye(2);

% Generate reference trajectory
x_mid_1 = [q_target(1)/4; 0.6;0;0;q_target(1)/problem_options.T;0;0;0];
x_mid_2 = [2*q_target(1)/4; 0.4;0;0;q_target(1)/problem_options.T;0;0;0];
x_mid_3 = [3*q_target(1)/4; 0.6;0;0;q_target(1)/problem_options.T;0;0;0];

x_target = [q_target;zeros(4,1)];
x_ref = interp1([0 0.25 0.5 0.75 1],[model.x0,x_mid_1,x_mid_2,x_mid_3,x_target]',linspace(0,1,settings.N_stages),'spline')'; %spline

model.lsq_x = {x, x_ref, Q}; % TODO also do trajectory
model.lsq_u = {u, u_ref, R}; % TODO also do trajectory
model.lsq_T = {x, x_target, Q_terminal};
%% terminal constraint and/or costs
%model.g_terminal = q(1:length(q_target))-q_target;

%% Solve OCP with NOSNOC

solver = NosnocSolver(model, settings);
if true
    x_guess = {};
    for ii = 1:settings.N_stages
        x_guess{ii} = x_ref(:,ii);
    end
    solver.set('x', x_guess');
end
[results,stats] = solver.solve();
%% plot
plot_results_hopping_robot

