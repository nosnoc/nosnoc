clear all
close all
import casadi.*

T = 1.0;
N_stages = 1;
N_finite_elements = 2;
%% Define (uncontrolled for now) projected system
model = nosnoc.model.Objects;
opts = nosnoc.Options;

%A1 = inv(diag([4, 0.1]));
A1 = inv(diag([1, 1]));
ball1 = nosnoc.objects.Ball(0.25,2);
ellipse2 = nosnoc.objects.Ellipse({A1});

ball1.x0 = [-4;-4];
ball1.x_dot = [1;1];

ellipse2.x0 = [0;0;0];
ellipse2.x_dot = [0;0;0];

model.addContact(ball1, ellipse2);

model.r = 0.01;

opts.T = T;
opts.N_stages = N_stages;
opts.N_finite_elements = N_finite_elements;
opts.n_s = 2;
opts.cross_comp_mode = CrossCompMode.FE_FE;
%opts.step_equilibration = 'mlcp';


%% Create problem

solver_opts = nosnoc.solver.Options();
%solver_opts.homotopy_steering_strategy = "ELL_INF";
%solver_opts.decreasing_s_elastic_upper_bound = true;
solver_opts.opts_casadi_nlp.ipopt.max_iter = 5000;
solver_opts.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
%solver_opts.opts_casadi_nlp.ipopt.acceptable_iter = 1;
%solver_opts.opts_casadi_nlp.ipopt.acceptable_tol = 1e-8;
solver_opts.warm_start_duals = true;
solver_opts.complementarity_tol = 1e-9;
%solver_opts.polishing_step = 1;
%solver_opts.N_homotopy = 1;
solver_opts.print_level = 5;

opts.preprocess();
model.verify_and_backfill(opts);

dcs = nosnoc.dcs.Objects(model);
dcs.generate_variables(opts);
dcs.generate_equations(opts);

dtp = nosnoc.discrete_time_problem.Objects(dcs, opts);

dtp.create_solver(solver_opts);

dtp.solve();
