clear all
close all
import casadi.*
import vdx.*

R = 3.5;
%% Define projected system
x1 = SX.sym('x1', 2);
x2 = SX.sym('x2', 2);
T = SX.sym('T');
x = [x1;x2];
x_target = [0;0;0;0];
data.x = [x;T];
data.lbx = [-inf;-inf;-inf;-inf;1e-1];
data.ubx = [inf;inf;inf;inf;inf];
x0 =[-25;-25;-15;-15];
data.x0 = [x0;1];
u1 = SX.sym('u1', 2);
u2 = SX.sym('u2', 2);
data.u = [u1;u2];
data.lbu = [-100/sqrt(2);-100/sqrt(2);-60/sqrt(2);-60/sqrt(2)];
data.ubu = [100/sqrt(2);100/sqrt(2);60/sqrt(2);60/sqrt(2)];
%data.lbu = [-100/sqrt(2);-100/sqrt(2);0;0];
%data.ubu = [100/sqrt(2);100/sqrt(2);0;0];
data.u0 = data.ubu;
data.c = [norm_2(x2-x1)-2*R];
data.f_x = T*[u1;u2;0];

% costs
data.f_q = 0;
data.f_q_T = 0.5*(x-x_target)'*(x-x_target) + 0.5*T^2;%0.5*(norm_2(x)^2);

data.T = 1;
data.N_stages = 25;
data.N_fe = 2;
data.n_s = 2;
data.irk_scheme = 'radau';

opts.step_eq = 'linear';
%opts.step_eq = 'direct_homotopy';
%opts.elastic_ell_inf = 1;

prob = InclusionProblem(data, opts);

prob.generate_constraints();

%% create solver
default_tol = 1e-12;

opts_casadi_nlp.ipopt.print_level = 1;
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

%% Do homotopy
prob.w.x(0,0,data.n_s).init = data.x0;
prob.w.x(0,0,data.n_s).lb = [x0;1e-1];
prob.w.x(0,0,data.n_s).ub = [x0;inf];
prob.w.lambda(0,0,data.n_s).init = 0;
prob.w.lambda(0,0,data.n_s).lb = 0;
prob.w.lambda(0,0,data.n_s).ub = 0;
prob.p.gamma_h(1).init = 1000;
homotopy(prob);
%% plot
fontsize = 12;
x_res = prob.w.x(0:data.N_stages,0:data.N_fe,data.n_s).res;
u_res = prob.w.u(1:data.N_stages).res;
u_rep = kron(u_res, ones(1,data.N_fe));
lambda_res = prob.w.lambda(1:data.N_stages, 0:data.N_fe, data.n_s).res;
h_res = prob.w.h(:,:).res;
t_res = [0,cumsum(h_res)]*x_res(end,1);
nabla_c_fun = casadi.Function('nabla_c', {data.x}, {data.c.jacobian(x)'});
nabla_c_res = full(nabla_c_fun(x_res));
c_res = full(prob.c_fun(x_res));
%plot_discs(h_res,x_res,[3.5,3.5], ["circle", "circle"])
v_res = u_rep + (1/x_res(end,1))*repmat(lambda_res,4,1).*nabla_c_res(:, 2:end);
figure('Position', [0,0, 400. 400])
ax = subplot(3,1,1);
plot(t_res(1:end-1), v_res([1,3],:), "LineWidth", 2)
xlabel("$t$")
ylabel("$v_x$")
ylim([40, 75])
xlim([0, x_res(end,1)])
grid on
ax.FontSize = fontsize;
ax.FontSize = fontsize;
ax = subplot(3,1,2);
plot(t_res(1:end-1), v_res([2,4],:), "LineWidth", 2)
xlabel("$t$")
ylabel("$v_y$")
ax.FontSize = fontsize;
ax.FontSize = fontsize;
ylim([40, 75])
xlim([0, x_res(end,1)])
grid on
ax = subplot(3,1,3);
plot(t_res, c_res, "LineWidth", 2)
xlabel("$t$")
ylabel("$c(x)$")
ax.FontSize = fontsize;
ax.FontSize = fontsize;
ylim([-0.1, 10])
xlim([0, x_res(end,1)])
grid on
