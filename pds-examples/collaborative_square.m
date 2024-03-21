clear all
close all
import casadi.*
import vdx.*

T = 10;
R = 0.5;
R_obj = 5;
%% Define projected system
x1 = SX.sym('x1', 2);
x2 = SX.sym('x2', 2);
x3 = SX.sym('x3', 3);
theta = x3(3);
R_matrix = [cos(theta) -sin(theta);...
           sin(theta) cos(theta)];
x = [x1;x2;x3];
data.x = x;
data.lbx = [-20;-20;-20;-20;-20;-20;-4*pi];
data.ubx = [20;20;20;20;20;20;4*pi];
data.x0 = [-10;0;0;10;0;0;0];
x_target = [-10;0;0;10;5;-10;0]; 
u1 = SX.sym('u1', 2);
u2 = SX.sym('u2', 2);
data.u = [u1;u2];
data.lbu = [-1;-1;-1;-1];
data.ubu = [1;1;1;1];
data.u0 = [0;0;0;0];
p = 4;
data.c = [sum((R_matrix*(x1-x3(1:2))).^p)-(R+R_obj)^p;sum((R_matrix*(x2-x3(1:2))).^p)-(R+R_obj)^p];
data.f_x = [u1;u2;0;0;0];

% costs
data.f_q = (x-x_target)'*diag([0,0,0,0,1,1,1])*(x-x_target);; %1e-1*norm_2(data.u)^2; %+ (x-x_target)'*diag([1e-6,1e-6,1e-6,1e-6,1,1,1e-6])*(x-x_target);
%data.f_q_T = (x-x_target)'*diag([1e-6,1e-6,1e-6,1e-6,1e3,1e3,1e3])*(x-x_target);
%data.f_q_T = (x-x_target)'*diag([0,0,0,0,1,1,1])*(x-x_target);
data.g_T = x(5:6)-x_target(5:6);

data.T = T;
data.N_stages = 20;
data.N_fe = 3;
data.n_s = 2;
data.irk_scheme = 'radau';

prob = InclusionProblem(data, struct);

prob.generate_constraints();

%% create solver
default_tol = 1e-7;

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
homotopy(prob,0.5);
%% plot
x_res = prob.w.x(0:data.N_stages,0:data.N_fe,data.n_s).res;
u_res = prob.w.u(1:data.N_stages).res;
h_res = prob.w.h(:).res;
t_res = [0,cumsum(h_res)];
plot_discs(h_res,x_res,[R,R,R_obj], ["circle", "circle", "square"])
