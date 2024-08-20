clear all;
clear all;
clc;
import casadi.*
%% init model and settings)
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
model = nosnoc.model.Cls();
%% Simulation settings
N_finite_elements = 4;
T_sim = 2;
N_sim = 21;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_finite_elements;
problem_options.N_sim = N_sim;
%%
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 2;
problem_options.use_fesd = 1;
problem_options.dcs_mode = 'CLS';
problem_options.friction_model = "Polyhedral";
problem_options.conic_model_switch_handling = "Abs";
problem_options.cross_comp_mode = 1;
problem_options.gamma_h = 1;
problem_options.fixed_eps_cls = 1;

solver_options.print_level = 3;
solver_options.N_homotopy = 10;
solver_options.homotopy_steering_strategy = 'ELL_INF';
solver_options.decreasing_s_elastic_upper_bound = true;
solver_options.use_previous_solution_as_initial_guess = 1;
%%
g = 10;
% Symbolic variables and bounds
q = SX.sym('q',3); 
v = SX.sym('v',3);
model.e = 0;
model.mu = 0.2;
model.dims.n_dim_contact = 3;
model.x = [q;v];
model.x0 = [0;0;1;2;1;0]; 
F_ext = [1;1]*0;
model.f_v = [F_ext;-g];
model.f_c = q(3);
model.J_tangent = [1 0;0 1; 0 0];
model.D_tangent = [1,-1,0,0;
                   0,0,1,-1;
                   0,0,0,0];
%% Call FESD Integrator
integrator = nosnoc.Integrator(model, problem_options, solver_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();
%%
qx = x_res(1,:);
qy = x_res(2,:);
qz = x_res(3,:);
vx = x_res(4,:);
vy = x_res(5,:);
vz = x_res(6,:);
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
%
figure
plot(t_grid,vx);
grid on
hold on
plot(t_grid,vy);
plot(t_grid ,vz);
xlabel('$t$','Interpreter','latex');
ylabel('$v$','Interpreter','latex');
legend({'$t_1^\top v$','$t_2^\top v$','$n^\top v$'},'Interpreter','latex','Location','best');
