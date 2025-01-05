clear all;
close all;
clc;
import casadi.*

%% 
J = 1; % no frictioinal impulse
J = 1/32; % frictional impulse apperas
above_ground = 0.1;
N_finite_elements = 5;
T_sim = 0.4;
N_sim = 41;
%% init nosnoc
problem_options = nosnoc.Options();
integrator_options = nosnoc.integrator.Options();
solver_options = integrator_options.fesd_solver_opts;
model = nosnoc.model.Cls();
%% settings
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 3;
problem_options.dcs_mode = 'Heaviside';
problem_options.pss_lift_step_functions = 1;
problem_options.time_freezing_Heaviside_lifting = 1;
problem_options.cross_comp_mode = 7;
problem_options.time_freezing = 1;
problem_options.impose_terminal_phyisical_time = 1;
problem_options.stagewise_clock_constraint = 0;
problem_options.a_n = 100;
problem_options.use_speed_of_time_variables = 1;
problem_options.local_speed_of_time_variable = 0;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_finite_elements;
problem_options.N_sim = N_sim;
problem_options.relax_terminal_physical_time = 'ELL_1';
problem_options.rho_terminal_physical_time = 1e6;
problem_options.s_sot_max = 100;

solver_options.homotopy_steering_strategy = 'ELL_INF';
solver_options.decreasing_s_elastic_upper_bound = true;
solver_options.opts_casadi_nlp.ipopt.max_iter = 1e3;
solver_options.print_level = 3;
solver_options.N_homotopy = 10;
integrator_opts.use_previous_solution_as_initial_guess = 0;

%%

model.e = 0;
model.mu = 1;
model.dims.n_dim_contact = 1;
%% the dynamics

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
m = 1; l = 1;
theta0 = pi/6;
g = 9.81;
M = diag([m,m,J]);
model.M = M;
% contact points of the rod
yc = qy-l/2*cos(qtheta);
xc = qx-l/2*sin(qtheta);
model.f_v = [0;-g;0];
model.f_c = yc;
model.J_tangent = xc.jacobian(q)';
% tangent
model.x0 = [0;l/2*cos(theta0)+above_ground;theta0 ;...
           -10;0;0];
%% Call FESD Integrator
integrator = nosnoc.Integrator(model, problem_options, integrator_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();
%%
qx = x_res(1,:);
qy = x_res(2,:);
qtheta = x_res(3,:);
vx = x_res(4,:);
vy = x_res(5,:);
omega = x_res(6,:);
t = x_res(7,:);
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
    pause(0.01)
    clf
end
%%
figure
plot(t,vx)
hold on
plot(t,vy)
plot(t,omega)
xlabel('$t$','Interpreter','latex')
ylabel('$v$','Interpreter','latex')
grid on
