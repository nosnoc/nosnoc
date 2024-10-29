clear all
close all
import casadi.*

T = 20.0;
N_stages = 30;
N_finite_elements = 4;
%% Define (uncontrolled for now) projected system
model = nosnoc.model.Objects;
opts = nosnoc.Options;
u = SX.sym('u', 2);
model.u = u;
model.lbu = [-1;-1];
model.ubu = [1;1];

A1 = inv(diag([4, 0.1]));
A2 = inv(diag([0.1, 4]));
ball1 = nosnoc.objects.Ball(0.25,2);
ellipse2 = nosnoc.objects.Ellipse({A1,A2});

ball1.x0 = [-3;-3];
ball1.x_dot = u;

ellipse2.x0 = [0;0;0];
ellipse2.x_dot = [0;0;0];

model.addContact(ball1, ellipse2);

model.r = 0.01;

opts.T = T;
opts.N_stages = N_stages;
opts.N_finite_elements = N_finite_elements;
opts.n_s = 4;
opts.cross_comp_mode = CrossCompMode.FE_FE;
%opts.step_equilibration = 'mlcp';
opts.gamma_h = 0.9;

model.f_q = u'*diag([1e-1,1e-1])*u;
ball1_target = [-3;-3];
ellipse_target = [0;0;2*pi];
x_target = vertcat(ball1_target, ellipse_target);
model.f_q_T = (model.x-x_target)'*diag([1e-1,1e-1,1e2,1e2,1e3])*(model.x-x_target);

%% Create problem

solver_opts = nosnoc.solver.Options();
%solver_opts.homotopy_steering_strategy = "ELL_INF";
%solver_opts.decreasing_s_elastic_upper_bound = true;
solver_opts.opts_casadi_nlp.ipopt.max_iter = 5000;
solver_opts.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
solver_opts.opts_casadi_nlp.ipopt.acceptable_iter = 1;
solver_opts.opts_casadi_nlp.ipopt.acceptable_tol = 1e-8;
solver_opts.warm_start_duals = true;
solver_opts.complementarity_tol = 1e-5;
solver_opts.polishing_step = 1;
solver_opts.print_level = 5;

model.verify_and_backfill(opts);

dcs = nosnoc.dcs.Objects(model);
dcs.generate_variables(opts);
dcs.generate_equations(opts);
