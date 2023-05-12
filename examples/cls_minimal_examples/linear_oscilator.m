clear all;
clear all;
clc;
import casadi.*
close all
%%
settings = NosnocOptions();
settings.irk_scheme = IRKSchemes.RADAU_IIA;
settings.print_level = 3;
settings.n_s = 2;
settings.N_homotopy = 5;
settings.cross_comp_mode = 1;
settings.dcs_mode = DcsMode.CLS;
settings.time_freezing = 0; %% we will need to exlude the coexistence of these two
settings.friction_model = "Conic";
% settings.friction_model = "Polyhedral";
settings.conic_model_switch_handling = "Plain";
settings.conic_model_switch_handling = "Abs";
% settings.conic_model_switch_handling = "Lp";
% settings.local_speed_of_time_variable = 0;
% settings.use_speed_of_time_variables = 0;
if settings.time_freezing
%     settings.dcs_mode = 'Step';
    settings.impose_terminal_phyisical_time = 1;
    settings.local_speed_of_time_variable = 1;
    settings.stagewise_clock_constraint = 0;
    settings.pss_lift_step_functions = 1;
    settings.sigma_0 = 1;
%     settings.step_equilibration = "heuristic_diff";
else
    settings.time_freezing = 0;
end
settings.mpcc_mode = "Scholtes_ineq";
settings.sigma_0  = 1e-2;

%%
g = 10;

% Symbolic variables and bounds
q = SX.sym('q',1);
v = SX.sym('v',1);
model.M = diag([1]);
model.x = [q;v];
model.e = 0;
model.mu = 0.0;
model.x0 = [0.2;-1];
k = 10; c = 0.3;
g = 0;
model.f_v = [-k*q-c*v-g+1];
model.f_c = q;
% model.J_tangent = [1; 0]; 
% model.D_tangent = [1,-1;0,0];

%% Simulation settings
N_FE = 5;
T_sim = 1;
N_sim = 1;
model.T_sim = T_sim;
model.N_FE = N_FE;
model.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 0;
%% Call nosnoc Integrator
[results,stats,model] = integrator_fesd(model,settings);
%% read and plot results
unfold_struct(results,'base');
qx = x_res(1,:);
vx = x_res(2,:);
if settings.time_freezing
t_opt = x_res(5,:);
else
    t_opt = results.t_grid;
end
figure
plot(t_opt,qx);
hold on
plot(t_opt,vx);
% axis equal
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');

% [qx(end),qy(end),t_opt(end)]

%%
if settings.time_freezing
alpha1 = alpha_res(1,:);
alpha2 = alpha_res(2,:);
alpha3 = alpha_res(3,:);
theta1 = alpha1+(1-alpha1).*(alpha2);
alpha_aux = (1-alpha1).*(1-alpha2);
theta2 = alpha_aux.*(1-alpha3);
theta3 = alpha_aux.*(alpha3);
figure;
subplot(131)
plot(t_grid,[theta1,nan])
xlabel('$\tau$','interpreter','latex');
ylabel(['$\theta_1$'],'interpreter','latex');
grid on
ylim([-0.1 1.1]);
subplot(132)
plot(t_grid,[theta2,nan])
xlabel('$\tau$','interpreter','latex');
ylabel(['$\theta_2$'],'interpreter','latex');
grid on
ylim([-0.1 1.1]);
subplot(133)
plot(t_grid,[theta3,nan])
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
stairs(s_sot_res)
grid on
xlabel('simulation step','interpreter','latex');
ylabel('$s$','interpreter','latex');
end
