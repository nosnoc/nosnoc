clear all;
close all
clc;
import casadi.*

%% init nosnoc
settings = NosnocOptions();  
model = NosnocModel();
%% settings
settings.irk_scheme = RKSchemes.RADAU_IIA;
settings.n_s = 1;
settings.mpcc_mode = 'elastic_ineq';
settings.opts_casadi_nlp.ipopt.max_iter = 1e3;
settings.print_level = 2;
settings.N_homotopy = 6;
settings.time_freezing = 1;
settings.pss_lift_step_functions = 1;
settings.stagewise_clock_constraint = 0;

%%
g = 10;
vertical_force = 0;
% Symbolic variables and bounds
q = SX.sym('q',2); 
v = SX.sym('v',2); 

model.x = [q;v]; 
model.e = 0;
model.mu_f = 0;
model.dims.n_dim_contact = 2;
model.a_n = g;
model.x0 = [0.8;0.5;-1.5;-1]; 
model.f_v = [0;-g];
model.f_c = [q(1);q(2)];
%% Simulation settings
N_FE = 3;
T_sim = 1.5;
N_sim = 40;
problem_options.T_sim = T_sim;
settings.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 0;
%% Call nosnoc Integrator
[results,stats,solver] = integrator_fesd(model,settings);
%% read and plot results
unfold_struct(results,'base');
qx = results.x(1,:); qy = results.x(2,:);
vx = results.x(3,:); vy = results.x(4,:);
t_opt = results.x(5,:);
figure
subplot(121)
plot(qx,qy,'LineWidth',2.5);
% axis equal
grid on
x_p = [-5 5 5 -5];
y_p = [-2 -2 0 0];
p = patch(x_p,y_p,'k');
p.FaceAlpha = 0.2; p.EdgeColor = 'none';
x_p = [-5 0 0 -5];
y_p = [0 0 5 5];
p = patch(x_p,y_p,'k');
p.FaceAlpha = 0.2; p.EdgeColor = 'none';
xlim([-1 2]);
ylim([-0.2 1]);
xlabel('$q_x$','interpreter','latex');
ylabel('$q_y$','interpreter','latex');
% axis equal
subplot(122)
plot(t_opt,vx);
hold on
plot(t_opt,vy);
legend({'$v_x$','$v_y$'},'interpreter','latex','Location','best');
% axis equal
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
