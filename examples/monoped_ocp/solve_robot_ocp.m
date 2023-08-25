function [varargout] = solve_robot_ocp(model,settings,scenario)
clc;
close all;
unfold_struct(scenario,'caller');
import casadi.*
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
model.u0 = u0;
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
use_unit_vectors = 1;
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
model.f_q = u'*u;
% box constraints
u_max = 200;
model.lbu = -u_max*ones(2,1);
model.ubu = u_max*ones(2,1);

% Sanity constraints
model.lbx = [-0.5;0;-pi;-pi;-inf;-inf;-inf;-inf];
model.ubx = [q_target(1)+0.5; 2.0;pi;pi;inf;inf;inf;inf];
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

% obstacle constraints
if hole_constraint
    temp = p_foot;
    p_foot_x = temp(1);
    p_foot_z = temp(2);
    a_vec = width_vec/2;
    b_vec = height_vec;
    for ii = 1:n_holes
        g_path = [g_path; ((p_foot_x-xc_vec(ii))/a_vec(ii))^2+((p_foot_z-zc_vec(ii))/b_vec(ii))^2-1];
    end
    g_path_lb = [g_path_lb;0*ones(n_holes ,1)];
    g_path_ub = [g_path_ub;inf*ones(n_holes ,1)];
end

if general_inequality_constraints
    if ~isempty(g_path)
        model.g_path = g_path;
        model.g_path_ub = g_path_ub;
        model.g_path_lb = g_path_lb;
    end
end

% least squares weight
Q = diag([1.0, 1.0, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8]);
Q_terminal = 100*diag([1.0, 1.0, 1e-8, 1e-8, 1.0, 1.0, 1.0, 1.0]);

% Generate reference trajectory
x_mid = [q_target(1)/2; 0.4;0;0;q_target(1)/problem_options.T;0;0;0];
x_target = [q_target;zeros(4,1)];
x_ref = interp1([0 0.5 1],[model.x0,x_mid,x_target]',linspace(0,1,settings.N_stages),'spline')' %spline

model.lsq_x = {x, x_ref, Q}; % TODO also do trajectory
model.lsq_T = {x, x_target, Q_terminal};
%% terminal constraint and/or costs
%model.g_terminal = q(1:length(q_target))-q_target;

%% Solve OCP with NOSNOC
if polishing_penalty_iteration
    solver = NosnocSolver(model, settings);
    [results,stats] = solver.solve();
%     w0 = [results.w_opt;stats.complementarity_stats(end)*10];
    w0 = [results.w_opt;1e-2];
    settings.s_elastic_max = 0.1;
%     settings.s_elastic_0 = 0.01;
    settings.mpcc_mode = 5;
    settings.N_homotopy = 1;
    settings.sigma_0 = 1e-6;
    solver = NosnocSolver(model, settings);
    [results,stats] = solver.solve();
else
    solver = NosnocSolver(model, settings);
    x_guess = {};
    for ii = 1:settings.N_stages
        x_guess{ii} = x_ref(:,ii);
    end
    solver.set('x', x_guess');
    [results,stats] = solver.solve();
end


%% Save statistics
fid = fopen('log_robot.txt','a');
fprintf(fid,[ '---------------------------------------\n']);
fprintf(fid,[ 'Scenario:' scenario.filename  '.\n']);
fprintf(fid,'Complementarity residual %2.2e \n',stats.complementarity_stats(end));
fprintf(fid,'CPU time %2.3f min. \n',stats.cpu_time_total/60);
fprintf(fid,[ '---------------------------------------\n']);
fclose(fid);
%%
results.stats = stats;
save(scenario.filename,'results')

%% plots and animation
plot_results_hopping_robot
%%
varargout{1} = results;
varargout{2} = stats;
varargout{3} = model;
varargout{4} = settings;

end


