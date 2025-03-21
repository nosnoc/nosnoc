clear all;
clear all;
clc;
import casadi.*

%% 
J = 1; % no frictioinal impulse
% J = 1/32; % frictional impulse apperas
above_ground = 0.2;

problem_options = nosnoc.Options();
integrator_options = nosnoc.integrator.Options();
solver_options = integrator_options.fesd_solver_opts;
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 2;
problem_options.dcs_mode = 'CLS';
problem_options.friction_model = "Polyhedral";
% problem_options.friction_model = "Conic";
% problem_options.conic_model_switch_handling = "Abs";
problem_options.pss_lift_step_functions = 1;


solver_options.opts_casadi_nlp.ipopt.max_iter = 3e3;
solver_options.print_level = 3;
solver_options.N_homotopy = 10;
solver_options.sigma_0 = 1e2;
solver_options.complementarity_tol = 1e-5;
solver_options.homotopy_steering_strategy = "ELL_INF";
% Some debug options
solver_options.print_details_if_infeasible = 0;
integrator_options.break_simulation_if_infeasible = 0;
%%
model = nosnoc.model.Cls();
model.e = 0;
model.mu = 0.2;

%% The CLS model
qx = SX.sym('qx',1);
qy = SX.sym('qy',1);
qtheta = SX.sym('qtheta',1);
vx = SX.sym('vx',1);
vy = SX.sym('vy',1);
omega = SX.sym('omega',1);
q = [qx;qy;qtheta];
v = [vx;vy;omega];
model.x = [q;v];
model.q = q;
model.v = v;
% constraint
m = 1; 
l = 1;
theta0 = pi/6;
g = 9.81;
M = diag([m,m,J]);
model.M = M;
% contact points of the rod
yc_left = qy-l/2*cos(qtheta);
xc_left = qx-l/2*sin(qtheta);
yc_right = qy+l/2*cos(qtheta);
xc_right = qx+l/2*sin(qtheta);
model.f_v = [0;-g;0];
model.f_c = [yc_left;yc_right];
model.J_tangent = [xc_left.jacobian(q)',xc_right.jacobian(q)'];
model.D_tangent = [xc_left.jacobian(q)',-xc_left.jacobian(q)'];

% tangent
% model.x0 = [0;l/2*cos(theta0)+above_ground;theta0 ;...
%            -10;0;0];
theta0 = pi/4;
model.x0 = [0;l/2*cos(theta0)+above_ground;theta0;...
           2;0;0];

%% Simulation settings
N_finite_elements = 2;
T_sim = 1;
N_sim = 10;

problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_finite_elements;
problem_options.N_sim = N_sim;
integrator_opts.use_previous_solution_as_initial_guess = 1;
%% Call FESD Integrator
integrator = nosnoc.Integrator(model, problem_options, integrator_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();

%% Plot results
qx = x_res(1,:);
qy = x_res(2,:);
qtheta = x_res(3,:);
vx = x_res(4,:);
vy = x_res(5,:);
omega = x_res(6,:);
if problem_options.time_freezing
    t = x_res(7,:);
else
    t = t_grid;
end

xc_res  = [];
yc_res  = [];
for ii = 1:length(qx)
    xc_res  = [xc_res, qx(ii)-l/2*sin(qtheta(ii))];
    yc_res  = [yc_res,qy(ii)-l/2*cos(qtheta(ii))];
end
%%
h = problem_options.h_k;
figure
for ii = 1:length(qx)
    plot([qx(ii)+l/2*sin(qtheta(ii)) xc_res(ii)],[qy(ii)+l/2*cos(qtheta(ii)) yc_res(ii)],'k','LineWidth',1.5)
    hold on
    yline(0,'r')
    xlabel('$q_x$','Interpreter','latex')
    ylabel('$q_y$','Interpreter','latex')
    axis equal
    ylim([-0.1 2])
    grid on
    pause(h)
    clf
end
%%
figure
subplot(311)
plot(t,vx)
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$v_x$','Interpreter','latex')
subplot(312)
plot(t,vy)
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$v_y$','Interpreter','latex')
subplot(313)
plot(t,omega)
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$\omega$','Interpreter','latex')

