clear all
close all
import casadi.*
import vdx.*

T = 5;
R = 1;
R_obj = 2;
%% Define projected system
x1 = SX.sym('x1', 2);
x2 = SX.sym('x2', 2);
x3 = SX.sym('x3', 2);
x = [x1;x2;x3];
x_target = [-10;0;10;0;0;10];
data.x = x;
data.lbx = [-inf;-inf;2.5;-inf;-inf;-inf];
data.ubx = [-2.5;inf;inf;inf;inf;inf];
data.x0 = [-10;10;5;-3;0;0];
u1 = SX.sym('u1', 2);
u2 = SX.sym('u2', 2);
data.u = [u1;u2];
data.lbu = [-10;-10;-10;-10];
data.ubu = [10;10;10;10];
data.u0 = [0;0;0;0];
data.c = [norm_2(x3-x1)-(R+R_obj);norm_2(x3-x2)-(R+R_obj)];
data.f_x = [u1;u2;0;0];

% costs
data.f_q = 1e-4*norm_2(data.u)^2;
data.f_q_T = (x-x_target)'*diag([1e0,1e0,1e0,1e0,1e3,1e3])*(x-x_target);

data.T = T;
data.N_stages = 50;
data.N_fe = 3;
data.n_s = 2;
data.irk_scheme = 'radau';

opts.step_eq = 'heuristic_mean';
opts.use_fesd = true;
%opts.elastic_ell_inf = 1;

prob = InclusionProblem(data, struct);

prob.generate_constraints();

%% create solver
default_tol = 1e-10;

%opts_casadi_nlp.ipopt.print_level = 1;
opts_casadi_nlp.print_time = 0;
opts_casadi_nlp.ipopt.sb = 'yes';
opts_casadi_nlp.verbose = false;
opts_casadi_nlp.ipopt.max_iter = 50000;
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

%% Do homotopy
prob.w.x(0,0,data.n_s).init = data.x0;
prob.w.x(0,0,data.n_s).lb = data.x0;
prob.w.x(0,0,data.n_s).ub = data.x0;
homotopy(prob);
%% plot
x_res = prob.w.x(0:data.N_stages,0:data.N_fe,data.n_s).res;
u_res = prob.w.u(1:data.N_stages).res;
h_res = prob.w.h(:).res;
t_res = [0,cumsum(h_res)];
fig = figure('Position', [10 10 1600 800]);
rectangle('Position',[-1.5 -25 3 50],'FaceColor',[1 0 0 0.2],'LineStyle', '--', 'EdgeColor' , 'red')
rectangle('Position',[-2 8 4 4], 'Curvature', 1, 'FaceColor',[0.8500 0.3250 0.0980,0.5], 'LineStyle', '--', 'EdgeColor' , [0.8500 0.3250 0.0980])
rectangle('Position',[-11 -1 2 2], 'Curvature', 1, 'FaceColor',[0 0.4470 0.7410,0.5], 'LineStyle', '--', 'EdgeColor' , [0 0.4470 0.7410])
rectangle('Position',[9 -1 2 2], 'Curvature', 1, 'FaceColor',[0 0.4470 0.7410,0.5], 'LineStyle', '--', 'EdgeColor' , [0 0.4470 0.7410])
plot_pass_discs(h_res,x_res,[R,R,R_obj], ["circle","circle","circle"], fig, 'coop');
