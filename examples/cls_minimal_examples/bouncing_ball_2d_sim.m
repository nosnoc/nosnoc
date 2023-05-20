clear all;
clc;
import casadi.*
close all
%%
settings = NosnocOptions();
settings.irk_scheme = IRKSchemes.RADAU_IIA;
settings.n_s = 2;
settings.print_level = 3;
settings.N_homotopy = 7;
settings.cross_comp_mode = 3;
settings.dcs_mode = DcsMode.CLS;
settings.no_initial_impacts = 1;
%settings.friction_model = "Polyhedral";
settings.friction_model = "Conic"; % "Conic"
settings.conic_model_switch_handling = "Abs";
% settings.mpcc_mode = MpccMode.elastic_ineq;
% settings.nlpsol = 'snopt';
% settings.psi_fun_type = CFunctionType.KANZOW_SCHWARTZ;
%%
g = 9.81;
% Symbolic variables and bounds
q = SX.sym('q',2);
v = SX.sym('v',2);
model = NosnocModel();
model.M = diag([1,1]);
model.x = [q;v];
model.e = 0;
model.mu = 0.2;
% impact
model.x0 = [0;1;4;0];
% Free flight
% model.x0 = [0;5.5;4;0];
% Friction slip slip
% model.x0 = [0;0.0;-0.6;0];
% friction slip stick
% model.x0 = [0;0.0;+0.6;0];
model.f_v = [3*1;-g];
model.f_c = q(2);
model.J_tangent = [1; 0];
model.D_tangent = [1,-1;0,0];
model.dims.n_dim_contact = 2; % TODO: REMOVE THIS IN time-freezing
%% Simulation settings
N_FE = 2;
T_sim = 1.7;
N_sim = 10;
model.T_sim = T_sim;
model.dims.N_finite_elements = N_FE;
model.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 1;
%% Call nosnoc Integrator
[results,stats,model] = integrator_fesd(model,settings);
%% read and plot results
unfold_struct(results,'base');
qx = results.x(1,:);
qy = results.x(2,:);
vx = results.x(3,:);
vy = results.x(4,:);
figure(1)
subplot(121)
plot(qx,qy);
axis equal
grid on
xlabel('$q_x$','interpreter','latex');
ylabel('$q_y$','interpreter','latex');
% axis equal
subplot(122)
plot(results.t_grid,vy);
hold on
plot(results.t_grid,vx);
% axis equal
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
%% forces
if N_sim == 1
figure
subplot(311)
plot(results.t_grid,results.x(3:4,:));
hold on
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
subplot(312)
plot(results.t_grid,[nan*ones(model.n_contacts,1),results.lambda_normal])
hold on
plot(results.t_grid,[nan*ones(model.n_tangents,1),results.lambda_tangent])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\lambda$','interpreter','latex');
subplot(313)
stem(results.t_grid,[results.Lambda_normal,nan*ones(model.n_contacts,1)])
hold on
stem(results.t_grid,[results.Lambda_tangent,nan*ones(model.n_tangents,1)]')
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda$','interpreter','latex');
xlim([0 T_sim])
end
