clear all
close all
import casadi.*
import vdx.*

T = 5.0;
N_sim = 100;
t_step = T/N_sim;
%% Define (uncontrolled for now) projected system
x = SX.sym('x', 2);
t = SX.sym('t', 1);
data.x = [x;t];
data.lbx = [-inf;-inf;-inf];
data.ubx = [inf;inf;inf];
data.x0 = [0;pi-2;0];
data.u = [];
data.lbu = [];
data.ubu = [];
data.u0 = [];
data.c = [-norm(x - [sin(t);cos(t)])^2+(1-0.5)^2];
data.f_x = [0; 0; 1];
data.f_q = 0;
data.f_q_T = 0;
data.partial_proj_matrix = diag([1 1 0]);

data.T = t_step;
data.N_stages = 1;
data.N_fe = 2;
data.n_s = 2;
data.irk_scheme = 'radau';

opts.step_eq = 'heuristic_mean';
opts.use_fesd = true;
opts.elastic_ell_inf = 0;
opts.time_dependent = true;

prob = InclusionProblem(data, opts);

prob.generate_constraints();

default_tol = 1e-12;

opts_casadi_nlp.ipopt.print_level = 2;
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
for step=1:N_sim
    prob.w.x(0,0,data.n_s).init = x_curr;
    prob.w.x(0,0,data.n_s).lb = x_curr;
    prob.w.x(0,0,data.n_s).ub = x_curr;
    success = homotopy(prob);
    if ~success
        disp(['Failure to converge at step=' num2str(step)])
    else
        disp(['step=' num2str(step)])
    end
    x_curr = prob.w.x(1,data.N_fe,data.n_s).res;
    x_sim = prob.w.x(1,:,data.n_s).res;
    x_sim_long = prob.w.x(1,:,:).res;
    %c_ind_sim = prob.w.c_ind(1,:,:).res';
    %c_ind_sim = [c_ind_sim{:}];
    lambda_sim = prob.w.lambda(1,:,:).res;
    if opts.use_fesd
        h_sim = prob.w.h(1,:).res;
    else
        h_sim = prob.p.T(1).init/(data.N_fe)*ones(data.N_fe,1)';
    end
    h_res = [h_res,h_sim];
    x_res = [x_res,x_sim];
    x_res_long = [x_res_long,x_sim_long];
    lambda_res = [lambda_res, lambda_sim];
    %c_ind = [c_ind,c_ind_sim];
    prob.w.init = prob.w.res;
end

c_res_long = full(c_fun(x_res_long));
c_res = full(c_fun(x_res));
t_res = [0,cumsum(h_res)];

fig = figure('Position', [10 10 1600 800]);
plot_moving_set(h_res,x_res,[0.5], ["circle"], fig, 'moving_set');

figure
hold on
plot(t_res, x_res(1,:))
plot(t_res, x_res(2,:))
hold off
if opts.use_fesd
    figure
    stairs(t_res(1:end-1),h_res)
end
figure
hold on
x1 = -1.7:0.001:0.7;
x2 = -(x1+0.5).^2 + 1;
plot(x_res(1,:), x_res(2,:))
hold off
