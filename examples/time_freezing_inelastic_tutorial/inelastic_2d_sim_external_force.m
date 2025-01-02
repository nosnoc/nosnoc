clear all;
clc;
import casadi.*
close all
%%
g = 10;
vertical_force = 0;
N_FE = 4;
T_sim = 1.5;
N_sim = 40;
u_sim = 1*ones(2,N_sim);

%% init nosnoc
problem_options = nosnoc.Options();
integrator_options = nosnoc.integrator.Options();
solver_options = integrator_options.fesd_solver_opts;
model = nosnoc.model.Cls();
%% settings
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 2;
problem_options.cross_comp_mode = 7;
problem_options.time_freezing = 1;
problem_options.pss_lift_step_functions = 0;
problem_options.impose_terminal_phyisical_time = 1;
problem_options.local_speed_of_time_variable = 1;
problem_options.stagewise_clock_constraint = 0;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;
problem_options.a_n = g;

solver_options.homotopy_steering_strategy = 'ELL_INF';
solver_options.decreasing_s_elastic_upper_bound = true;
solver_options.opts_casadi_nlp.ipopt.max_iter = 1e3;
solver_options.print_level = 3;
solver_options.N_homotopy = 6;
integrator_opts.use_previous_solution_as_initial_guess = 0;
%%
% Symbolic variables and bounds
q = SX.sym('q',2); v = SX.sym('v',2); 
u = SX.sym('u',2);

model.x = [q;v];
model.u = u;
model.e = 0;
model.mu = 0.3;
model.dims.n_dim_contact = 1;
model.x0 = [0;1;3;0]; 
model.f_v = [0;-g+vertical_force*g*q(1)]+u;
model.f_c = q(2);
model.J_tangent = [1; 0];

%% Call nosnoc Integrator
integrator = nosnoc.Integrator(model, problem_options, integrator_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate("u", u_sim);

%% read and plot results
qx = x_res(1,:);
qy = x_res(2,:);
vx = x_res(3,:);
vy = x_res(4,:);
t_opt = x_res(5,:);
figure
subplot(121)
plot(qx,qy);
axis equal
grid on
xlabel('$q_x$','interpreter','latex');
ylabel('$q_y$','interpreter','latex');
% axis equal
subplot(122)
plot(t_opt,vy);
hold on
plot(t_opt,vx);
% axis equal
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');

% [qx(end),qy(end),t_opt(end)]

%%
alpha = integrator.get('alpha');
alpha1 = alpha(1,:);
alpha2 = alpha(2,:);
alpha3 = alpha(3,:);
theta1 = alpha1+(1-alpha1).*(alpha2);
alpha_aux = (1-alpha1).*(1-alpha2);
theta2 = alpha_aux.*(1-alpha3);
theta3 = alpha_aux.*(alpha3);
figure;
subplot(131)
plot(t_grid,theta1)
xlabel('$\tau$','interpreter','latex');
ylabel(['$\theta_1$'],'interpreter','latex');
grid on
ylim([-0.1 1.1]);
subplot(132)
plot(t_grid,theta2)
xlabel('$\tau$','interpreter','latex');
ylabel(['$\theta_2$'],'interpreter','latex');
grid on
ylim([-0.1 1.1]);
subplot(133)
plot(t_grid,theta3)
xlabel('$\tau$','interpreter','latex');
ylabel(['$\theta_3$'],'interpreter','latex');
grid on
ylim([-0.1 1.1]);
%% speed of time
figure
subplot(121)
plot(t_grid,t_opt)
hold on
plot(t_grid,t_grid,'k--')
grid on
xlabel('$\tau$','interpreter','latex');
ylabel('$t$','interpreter','latex');
subplot(122)
s_sot_res = integrator.get('sot');
stairs(s_sot_res)
grid on
xlabel('simulation step','interpreter','latex');
ylabel('$s$','interpreter','latex');



%% complementarity residuals
if 0 
lambda0 = lambda_0_res;
lambda1 = lambda_1_res;
alpha = alpha;
comp1 = alpha.*lambda0;
comp2 = lambda1.*(ones(size(alpha))-alpha);
figure
subplot(121)
plot(comp1')
subplot(122)
plot(comp2')
end
