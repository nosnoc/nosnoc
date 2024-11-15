clear all;
close all
clc;
import casadi.*

%% init nosnoc
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
model = nosnoc.model.Cls();
%%
g = 10;
N_FE = 3;
T_sim = 1.5;
N_sim = 40;
%% settings
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 1;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;
problem_options.time_freezing = 1;
problem_options.pss_lift_step_functions = 1;
problem_options.stagewise_clock_constraint = 0;
problem_options.a_n = g;


solver_options.homotopy_steering_strategy = 'ELL_INF';
solver_options.opts_casadi_nlp.ipopt.max_iter = 1e3;
solver_options.print_level = 2;
solver_options.N_homotopy = 6;
solver_options.use_previous_solution_as_initial_guess = 0;

%%
vertical_force = 0;
% Symbolic variables and bounds
q = SX.sym('q',2); 
v = SX.sym('v',2); 

model.x = [q;v]; 
model.e = 0;
model.mu = 0;
model.dims.n_dim_contact = 2;
model.x0 = [0.8;0.5;-1.5;-1]; 
model.f_v = [0;-g];
model.f_c = [q(1);q(2)];
%% Call nosnoc Integrator
integrator = nosnoc.integrator.FESD(model, problem_options, solver_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();
%% read and plot results
qx = x_res(1,:); qy = x_res(2,:);
vx = x_res(3,:); vy = x_res(4,:);
t_opt = x_res(5,:);
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
