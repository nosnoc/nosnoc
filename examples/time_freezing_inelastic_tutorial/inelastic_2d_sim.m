%clear all;
clc;
import casadi.*
close all
%%
g = 10;
vertical_force = 0;
%% init nosnoc
problem_options = nosnoc.Options();
integrator_options = nosnoc.integrator.Options();
solver_options = integrator_options.fesd_solver_opts;
model = nosnoc.model.Cls();
%%
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
solver_options.print_level = 3;
problem_options.cross_comp_mode = 7;
%problem_options.lift_complementarities = 1;
problem_options.time_freezing = 1;
problem_options.impose_terminal_phyisical_time = 1;
problem_options.local_speed_of_time_variable = 1;
problem_options.stagewise_clock_constraint = 0;
solver_options.complementarity_tol = 1e-10;
solver_options.sigma_N = 1e-10;
solver_options.N_homotopy = 16;
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
problem_options.a_n = 10;

problem_options.pss_lift_step_functions = 0;
%%
% Symbolic variables and bounds
q = SX.sym('q',2); 
v = SX.sym('v',2);

model.x = [q;v]; 
model.e = 0;
model.mu = 0.0;
model.x0 = [0;0.06;3;0]; 
model.f_v = [0;-g+vertical_force*g*q(1)];
model.f_c = q(2);
model.J_tangent = [1; 0];
model.dims.n_dim_contact = 2;

%% Simulation settings
N_FE = 3;
T_sim = 0.2;
N_sim = 1;
problem_options.T_sim = T_sim;
problem_options.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;
problem_options.n_s = 2;
integrator_opts.use_previous_solution_as_initial_guess = 0;
%% Call nosnoc Integrator
integrator = nosnoc.Integrator(model, problem_options, integrator_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();
results.x = x_res;
results.alpha = integrator.get("alpha");
results.t_grid = t_grid;
results.s_sot = integrator.get("sot");
%% read and plot results
plot_particle(results, model.mu ~= 0)

figure
subplot(411)
plot(integrator.get_full("lambda_n")')
ylabel('$\lambda_n$');
subplot(412)
plot(integrator.get_full("alpha")')
ylabel('$\alpha$');
subplot(413)
plot(integrator.get_full("lambda_p")')
ylabel('$\lambda_p$');
subplot(414)
plot(1-integrator.get_full("alpha")')
ylabel('$1-\alpha$');

figure
z = integrator.get_full("z");
subplot(311)
plot(z(1,:))
ylabel("$\theta_1$")
subplot(312)
plot(z(2,:))
ylabel("$\theta_2$")
subplot(313)
plot(z(3,:))
ylabel("$\beta$")
