clear all
close all
import casadi.*

T_sim = 5.0;
N_sim = 5;
N_finite_elements = 2;
%% Define (uncontrolled for now) projected system
model = nosnoc.model.PDSObjects;
problem_options = nosnoc.Options;
integrator_options = nosnoc.integrator.Options();
solver_options = integrator_options.fesd_solver_opts;

%A1 = inv(diag([4, 0.1]));
A1 = inv(diag([1, 1]));
ball1 = nosnoc.objects.Ball(0.25,2);
ellipse2 = nosnoc.objects.Ellipse({A1});

ball1.x0 = [-4;-4];
ball1.x_dot = [1;1];

ellipse2.x0 = [0;0;0];
ellipse2.x_dot = [0;0;0];

model.addContact(ball1, ellipse2, 0.1);

problem_options.N_finite_elements = N_finite_elements;
problem_options.n_s = 2;
problem_options.cross_comp_mode = CrossCompMode.FE_FE;

problem_options.N_sim = N_sim;
problem_options.T_sim = T_sim;
%opts.step_equilibration = 'mlcp';


%% Create problem
solver_options.homotopy_steering_strategy = "ELL_INF";
solver_options.decreasing_s_elastic_upper_bound = true;
solver_options.opts_casadi_nlp.ipopt.max_iter = 5000;
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
%solver_options.opts_casadi_nlp.ipopt.acceptable_iter = 1;
%solver_options.opts_casadi_nlp.ipopt.acceptable_tol = 1e-8;
solver_options.warm_start_duals = true;
solver_options.complementarity_tol = 1e-9;
%solver_options.polishing_step = 1;
%solver_options.N_homotopy = 1;
solver_options.print_level = 3;

integrator = nosnoc.Integrator(model, problem_options, integrator_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();

