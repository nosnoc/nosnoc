clear all;
clc;
import casadi.*
close all
%% init nosnoc
settings = NosnocOptions();
model = nosnoc.model.Cls();
%% settings
settings.rk_scheme = RKSchemes.RADAU_IIA;
settings.n_s = 1;
settings.print_level = 2;
settings.N_homotopy = 6;
settings.cross_comp_mode = 1;
settings.time_freezing = 1;
settings.pss_lift_step_functions = 0;
settings.impose_terminal_phyisical_time = 1;
settings.local_speed_of_time_variable = 1;
settings.stagewise_clock_constraint = 0;
%%
g = 10;
vertical_force = 0;
% Symbolic variables and bounds
q = SX.sym('q',2); v = SX.sym('v',2); 
u = SX.sym('u',2);

model.x = [q;v];
model.u = u;
model.e = 0;
model.mu_f = 0.3;
model.dims.n_dim_contact = 2;
model.a_n = g;
model.x0 = [0;1;3;0]; 
model.f_v = [0;-g+vertical_force*g*q(1)]+u;
model.f_c = q(2);
model.J_tangent = [1; 0];

%% Simulation settings
N_FE = 3;
T_sim = 1.5;
N_sim = 40;
u_sim = 1*ones(2,N_sim);
problem_options.T_sim = T_sim;
settings.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 0;
%% Call nosnoc Integrator
[results,stats] = integrator_fesd(model,settings,u_sim);
unfold_struct(results,'base');

%% read and plot results
qx = results.x(1,:);
qy = results.x(2,:);
vx = results.x(3,:);
vy = results.x(4,:);
t_opt = results.x(5,:);
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
alpha1 = results.alpha(1,:);
alpha2 = results.alpha(2,:);
alpha3 = results.alpha(3,:);
theta1 = alpha1+(1-alpha1).*(alpha2);
alpha_aux = (1-alpha1).*(1-alpha2);
theta2 = alpha_aux.*(1-alpha3);
theta3 = alpha_aux.*(alpha3);
figure;
subplot(131)
plot(results.t_grid,[theta1,nan])
xlabel('$\tau$','interpreter','latex');
ylabel(['$\theta_1$'],'interpreter','latex');
grid on
ylim([-0.1 1.1]);
subplot(132)
plot(results.t_grid,[theta2,nan])
xlabel('$\tau$','interpreter','latex');
ylabel(['$\theta_2$'],'interpreter','latex');
grid on
ylim([-0.1 1.1]);
subplot(133)
plot(results.t_grid,[theta3,nan])
xlabel('$\tau$','interpreter','latex');
ylabel(['$\theta_3$'],'interpreter','latex');
grid on
ylim([-0.1 1.1]);
%% speed of time
figure
subplot(121)
plot(results.t_grid,t_opt)
hold on
plot(results.t_grid,results.t_grid,'k--')
grid on
xlabel('$\tau$','interpreter','latex');
ylabel('$t$','interpreter','latex');
subplot(122)
stairs(s_sot_res)
grid on
xlabel('simulation step','interpreter','latex');
ylabel('$s$','interpreter','latex');



%% complementarity residuals
if 0 
lambda0 = lambda_0_res;
lambda1 = lambda_1_res;
alpha = results.alpha;
comp1 = alpha.*lambda0;
comp2 = lambda1.*(ones(size(alpha))-alpha);
figure
subplot(121)
plot(comp1')
subplot(122)
plot(comp2')
end
