clear all
close all
import casadi.*

T = 20.0;
N_stages = 40;
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
opts.n_s = 2;
opts.cross_comp_mode = CrossCompMode.FE_FE;
%opts.step_equilibration = 'mlcp';
opts.gamma_h = 0.9;

model.f_q = u'*diag([1e-1,1e-1])*u;
ball1_target = [-3;-3];
ellipse_target = [1;1;2*pi];
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
solver_opts.complementarity_tol = 1e-9;
%solver_opts.polishing_step = 1;
solver_opts.print_level = 5;

opts.preprocess();
model.verify_and_backfill(opts);

dcs = nosnoc.dcs.Objects(model);
dcs.generate_variables(opts);
dcs.generate_equations(opts);

dtp = nosnoc.discrete_time_problem.Objects(dcs, opts);

dtp.create_solver(solver_opts);

dtp.solve();

x_res = dtp.w.x(0:opts.N_stages,0:N_finite_elements,opts.n_s).res;
u_res = dtp.w.u(1:opts.N_stages).res;
pd_res = dtp.w.p_d(0:opts.N_stages,0:N_finite_elements,opts.n_s).res;
h_res = dtp.w.h.res;

pgon1 = ball1.to_polygon();
pgon2 = ellipse2.to_polygon(model.r);

facecolor1 = [0 0.4470 0.7410];
linecolor1 = facecolor1*0.7;
facecolor2 = [0.8500 0.3250 0.0980];
linecolor2 = facecolor2*0.7;

fig = figure('Position', [10 10 1600 800]);

indices = {1:2, 3:5};

plot_pds_sdf_example(h_res, x_res, pd_res, indices, [pgon1,pgon2], {facecolor1,facecolor2}, {linecolor1,linecolor2,}, fig, 'ellipse_spinning_friction');

