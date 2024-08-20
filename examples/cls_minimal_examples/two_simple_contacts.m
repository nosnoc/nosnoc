clear all;
clear all;
clc;
import casadi.*
close all
%%
settings = NosnocOptions();
model = nosnoc.model.Cls();
%%
settings.rk_scheme = RKSchemes.GAUSS_LEGENDRE;
settings.n_s = 2;
settings.print_level = 3;
settings.N_homotopy = 20;
settings.mpcc_mode = 'elastic_ineq';
settings.decreasing_s_elastic_upper_bound = 1;
settings.cross_comp_mode = 1;
settings.dcs_mode = DcsMode.CLS;
% settings.friction_model = "Polyhedral";
settings.friction_model = "Conic"; % "Conic"
settings.conic_model_switch_handling = "Abs";
settings.no_initial_impacts = 1;
%%
g = 9.81;
q = SX.sym('q',2);
v = SX.sym('v',2);
model.M = diag([1,1]);
model.x = [q;v];
model.e = [1;0];
model.mu_f = [0.4;0.7]*0;
model.x0 = [0.2;0.3;-2;0];
model.f_v = [3*1;-g];
model.f_c = [q(1);q(2)];
model.J_tangent = [0 1;1 0];
model.D_tangent = [model.J_tangent -model.J_tangent] ;
%% Simulation settings
N_FE = 5;
T_sim = 1;
N_sim = 1;
problem_options.T_sim = T_sim;

settings.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 1;
%% Call nosnoc Integrator
[results,stats,solver] = integrator_fesd(model,settings);
%% read and plot results
unfold_struct(results,'base');
qx = results.x(1,:);
qy = results.x(2,:);
vx = results.x(3,:);
vy = results.x(4,:);
t_grid = results.t_grid;
figure(1)
subplot(121)
plot(qx,qy);
axis equal
grid on
hold on
xline(0,'k')
yline(0,'k')
xlim([-0.2 1])
ylim([-0.2 1])
xlabel('$q_x$','interpreter','latex');
ylabel('$q_y$','interpreter','latex');
% axis equal
subplot(122)
plot(t_grid,vy);
hold on
plot(t_grid,vx);
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
%% forces
figure
subplot(311)
plot(t_grid,results.x(3:4,:));
hold on
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
subplot(312)
plot(t_grid,[nan*ones(model.dims.n_contacts,1),results.lambda_normal])
hold on
plot(t_grid,[nan*ones(model.dims.n_tangents,1),results.lambda_tangent])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\lambda$','interpreter','latex');
subplot(313)
stem(t_grid,[results.Lambda_normal,nan*ones(model.dims.n_contacts,1)]')
hold on
stem(t_grid,[results.Lambda_tangent,nan*ones(model.dims.n_tangents,1)]')
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda$','interpreter','latex');
xlim([0 T_sim])

