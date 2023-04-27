clear all;
clear all;
clc;
import casadi.*
close all
%%
settings = NosnocOptions();
settings.irk_scheme = IRKSchemes.RADAU_IIA;
settings.n_s = 2;
% settings.psi_fun_type = CFunctionType.STEFFENSON_ULBRICH;
settings.print_level = 2;
settings.N_homotopy = 6;
settings.cross_comp_mode = 6;
settings.dcs_mode = DcsMode.CLS;
settings.time_freezing = 0; %% we will need to exlude the coexistence of these two
settings.friction_model = "Conic"; % "Polyhedral"
settings.conic_model_switch_handling = "Abs";  % Plain % Lp
settings.local_speed_of_time_variable = 0;
settings.use_speed_of_time_variables = 0;
% if settings.time_freezing && settings.dcs_mode~=DcsMode.CLS
% settings.impose_terminal_phyisical_time = 1;
% settings.local_speed_of_time_variable = 1;
% settings.stagewise_clock_constraint = 0;
% settings.pss_lift_step_functions = 0;
% else
%     settings.time_freezing = 1;
% end
%%
g = 10;
% Symbolic variables and bounds
q = SX.sym('q',1);
v = SX.sym('v',1);
model.M = diag([1]);
model.x = [q;v];
model.e = 1;
model.mu = 0;
model.a_n = 20;
model.x0 = [0.2;0];
model.f_v = [-g];
model.f_c = q;
model.J_tangent = [1; 0]; 
model.D_tangent = [1,-1;0,0];
%% Simulation setings
N_FE = 2;
T_sim = 0.8*10;
N_sim = 20;
model.T_sim = T_sim;
model.N_FE = N_FE;
model.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 1;
%% Call nosnoc Integrator
[results,stats,model,settings,solver] = integrator_fesd(model,settings);
%% read and plot results
unfold_struct(results,'base');
qx = x_res(1,:);
vx = x_res(2,:);
%t_opt = x_res(5,:);
t_opt = t_grid;
figure
subplot(121)
plot(t_opt,qx);
axis equal
grid on
xlabel('$q_x$','interpreter','latex');
ylabel('$t$','interpreter','latex');
% axis equal
subplot(122)
plot(t_opt,vx);
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');

% [qx(end),qy(end),t_opt(end)]

