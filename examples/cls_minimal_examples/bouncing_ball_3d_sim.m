clear all;
clear all;
clc;
import casadi.*
%%
[settings] = NosnocOptions();  
settings.irk_scheme = IRKSchemes.RADAU_IIA;
settings.n_s = 2;
settings.mpcc_mode = MpccMode.Scholtes_ineq;
settings.print_level = 3;
settings.N_homotopy = 10;
settings.use_fesd = 1;
%settings.time_freezing = 1;
settings.dcs_mode = 'CLS';
settings.friction_model = "Polyhedral";
%settings.friction_model = "Conic";
settings.conic_model_switch_handling = "Abs";
settings.pss_lift_step_functions= 0;
settings.impose_terminal_phyisical_time  = 1;
settings.stagewise_clock_constraint = 0;
settings.nonsmooth_switching_fun = 0;
settings.pss_lift_step_functions = 0;
settings.cross_comp_mode = 1;
%%
g = 10;
% Symbolic variables and bounds
q = SX.sym('q',3); 
v = SX.sym('v',3);
model = NosnocModel();
model.e = 0;
model.mu = 0.2;
model.dims.n_dim_contact = 3;
model.x = [q;v]; 
model.a_n = g;
model.x0 = [0;0;1;2;1;0]; 
F_ext = [1;1]*0;
% norm(F_ext)
model.f_v = [F_ext;-g];
model.f_c = q(3);
model.J_tangent = [1 0;0 1; 0 0];
model.D_tangent = [1,-1,0,0;
                   0,0,1,-1;
                   0,0,0,0];
%% Simulation settings
N_finite_elements = 10;
T_sim = 2;
N_sim = 1;
model.T_sim = T_sim;
settings.N_finite_elements = N_finite_elements;
model.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 0;
%% Call FESD Integrator
[results,stats,solver] = integrator_fesd(model,settings);
%%
qx = results.x(1,:);
qy = results.x(2,:);
qz = results.x(3,:);
vx = results.x(4,:);
vy = results.x(5,:);
vz = results.x(6,:);
t_grid = results.t_grid;
figure
plot3(qx,qy,qz);
axis equal
xlim([-0.1 2])
ylim([-0.1 2])
zlim([-0.1 1])
grid on
xlabel('$q_x$','Interpreter','latex');
ylabel('$q_y$','Interpreter','latex');
zlabel('$q_z$','Interpreter','latex');
%%
figure
plot(t_grid,vx,'LineWidth',2);
grid on
hold on
plot(t_grid,vy);
plot(t_grid ,vz);
xlabel('$t$','Interpreter','latex');
ylabel('$v$','Interpreter','latex');
legend({'$t_1^\top v$','$t_2^\top v$','$n^\top v$'},'Interpreter','latex','Location','best');
