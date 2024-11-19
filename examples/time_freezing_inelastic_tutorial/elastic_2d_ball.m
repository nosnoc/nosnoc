clear all;
clear all;
clc;
import casadi.*
close all
%%
problem_options = nosnoc.Options();
integrator_options = nosnoc.integrator.Options();
solver_options = integrator_options.fesd_solver_opts;
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.dcs_mode = DcsMode.Stewart;
problem_options.n_s = 3;
problem_options.cross_comp_mode = 'FE_FE';
problem_options.use_fesd = 1;
problem_options.time_freezing = 1;
problem_options.pss_lift_step_functions = 0;


problem_options.impose_terminal_phyisical_time = 1;
problem_options.local_speed_of_time_variable = 0;
problem_options.stagewise_clock_constraint = 0;
problem_options.k_aux = 100;
problem_options.relax_terminal_numerical_time = 'ELL_2';
problem_options.relax_terminal_physical_time = 'ELL_2';
problem_options.rho_terminal_physical_time = 1e6;


solver_options.complementarity_tol = 1e-5;
solver_options.opts_casadi_nlp.ipopt.max_iter = 1e4;
solver_options.homotopy_steering_strategy = 'ELL_INF';
solver_options.decreasing_s_elastic_upper_bound = true;
%%
g = 10;
% Symbolic variables and bounds
q = SX.sym('q',1); 
v = SX.sym('v',1);
model = nosnoc.model.Cls();
model.x = [q;v]; 
model.e = 0.9;
model.x0 = [1;0]; 
model.f_v = -g;
model.f_c = q;
model.dims.n_dim_contact = 1;
%% Simulation settings
N_FE = 5;
T_sim = 2;
N_sim = 61;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;
solver_options.use_previous_solution_as_initial_guess = 1;
solver_options.print_level = 3;
%% Call nosnoc Integrator
integrator = nosnoc.Integrator(model, problem_options, integrator_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();
%% read and plot results
qx = x_res(1,:);
vx = x_res(2,:);
t_opt = x_res(3,:);
t_grid = integrator.get_time_grid();
s_sot = integrator.get('sot');
figure
subplot(121)
plot(t_opt,qx);
axis equal
grid on
xlabel('$t$','interpreter','latex');
ylabel('$q$','interpreter','latex');
% axis equal
subplot(122)
plot(t_opt,vx);
hold on
% axis equal
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');

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
stairs(s_sot)
grid on
xlabel('simulation step','interpreter','latex');
ylabel('$s$','interpreter','latex');

%% complementarity residuals
