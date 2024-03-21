clear all
close all
import casadi.*
import vdx.*

T = 0.5;
%% Define (uncontrolled for now) projected system
x = SX.sym('x', 2);
data.x = x;
data.lbx = [-inf;-inf];
data.ubx = [inf;inf];
data.x0 = [0.23;-1];
%data.x0 = [0.8;-.8];
data.u = [];
data.lbu = [];
data.ubu = [];
data.u0 = [];
data.c = [x(2)+1];
data.f_x = [x(2); -x(1)];
data.f_q = 0;
data.f_q_T = 0;

data.T = T;
data.N_stages = 1;
data.N_fe = 3;
data.n_s = 2;
data.irk_scheme = 'radau';

opts.step_eq = 'linear';
opts.use_fesd = true;
%opts.elastic_ell_inf = 1;

prob = InclusionProblem(data, opts);

prob.generate_constraints();

default_tol = 1e-12;

opts_casadi_nlp.ipopt.print_level = 0;
opts_casadi_nlp.print_time = 0;
opts_casadi_nlp.ipopt.sb = 'yes';
opts_casadi_nlp.verbose = false;
opts_casadi_nlp.ipopt.max_iter = 10000;
opts_casadi_nlp.ipopt.bound_relax_factor = 0;
%opts_casadi_nlp.ipopt.bound_relax_factor = 1e-8;
%opts_casadi_nlp.ipopt.honor_original_bounds = 'yes';
opts_casadi_nlp.ipopt.tol = default_tol;
opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
opts_casadi_nlp.ipopt.compl_inf_tol = default_tol;
opts_casadi_nlp.ipopt.acceptable_tol = 1e-6;
opts_casadi_nlp.ipopt.mu_strategy = 'adaptive';
opts_casadi_nlp.ipopt.mu_oracle = 'quality-function';
opts_casadi_nlp.ipopt.warm_start_init_point = 'yes';
opts_casadi_nlp.ipopt.warm_start_entire_iterate = 'yes';
opts_casadi_nlp.ipopt.linear_solver = 'ma27';
prob.create_solver(opts_casadi_nlp);

c_fun = Function('c_fun', {data.x}, {data.c});

x_res = data.x0;
x_res_long = data.x0;
c_ind = [0;0];
h_res = [];
lambda_res = [];
x_curr = data.x0;
prob.w.x(0,0,data.n_s).init = x_curr;
prob.w.x(0,0,data.n_s).lb = x_curr;
prob.w.x(0,0,data.n_s).ub = x_curr;
status = homotopy(prob);
