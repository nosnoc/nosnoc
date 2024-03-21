clear all
close all
import casadi.*
import vdx.*

%% Define (uncontrolled for now) projected system
x = SX.sym('x', 2);
c = [x(2)+0.1;-x(2) - (x(1)+0.5)^2 + 1];
f = [x(2); -x(1)];

%% Define the needed functions
n_x = length(x);
n_c = length(c);

lambda = SX.sym('lambda', n_c);

nabla_c = c.jacobian(x)';

f_x = f + nabla_c*lambda;

f_x_fun = Function('f_x_fun', {x,lambda}, {f_x});
c_fun = Function('c_fun', {x}, {c});

%% FESD problem data
N_sim = 1000;
N_stages = 1;
N_fe = 3;
n_s = 2;

T_sim = 10.0;
T = T_sim/N_sim;
t_stage = T/N_stages;
h0 = t_stage/N_fe;
irk_scheme = 'radau';
[B, C, D, tau_root] = generate_butcher_tableu_integral(n_s, irk_scheme);

prob = Problem();

%% Create variables for problem
prob.w.x(0,0,0) = {{['x_0'], n_x}};
for ii=1:N_stages
    for jj=1:N_fe
        prob.w.h(ii,jj) = {{['h_' num2str(ii) '_' num2str(jj)], 1}, 0, 2*h0, h0};
        for kk=1:n_s
            prob.w.x(ii,jj,kk) = {{['x_' num2str(ii) '_' num2str(jj) '_' num2str(kk)], n_x}};
            prob.w.lambda(ii,jj,kk) = {{['lambda_' num2str(ii) '_' num2str(jj) '_' num2str(kk)], n_c},0,inf};
        end
    end
end

prob.p.sigma(1) = {SX.sym('sigma'),0,inf,0};
prob.p.gamma_h(1) = {SX.sym('gamma_h'),0,inf,0.1};

%% Dyanamics
x_prev = prob.w.x(0,0,0);
for ii=1:N_stages
    for jj=1:N_fe
        for kk=1:n_s
            x_ijk = prob.w.x(ii,jj,kk);
            lambda_ijk = prob.w.lambda(ii,jj,kk);
            fj = f_x_fun(x_ijk,lambda_ijk);
            xk = C(1, kk+1) * x_prev;
            for rr=1:n_s
                x_ijr = prob.w.x(ii,jj,rr);
                xk = xk + C(rr+1, kk+1) * x_ijr;
            end
            h = prob.w.h(ii,jj);
            prob.g.dynamics(ii,jj,kk) = {h * fj - xk};
            % also add non-negativity constraint on c
            prob.g.c_nonnegative(ii,jj,kk) = {c_fun(x_ijk), 0, inf};
        end
        x_prev = prob.w.x(ii,jj,n_s);
    end
end

%% Sum of hs
for ii=1:N_stages
    sum_h = 0;
    for jj=1:N_fe
        sum_h = sum_h + prob.w.h(ii,jj);
    end
    prob.g.sum_h(ii) = {t_stage-sum_h};
end

%% Cross Complementarity
% In this case using FE-FE only because its easy to implement :)

x_prev = prob.w.x(0,0,0);
G = [];
H = [];
for ii=1:N_stages
    for jj=1:N_fe
        Gij = c_fun(x_prev);
        Hij = 0;
        for kk=1:n_s
            x_ijk = prob.w.x(ii,jj,kk);
            lambda_ijk = prob.w.lambda(ii,jj,kk);
            Gij = Gij + c_fun(x_ijk);
            Hij = Hij + lambda_ijk;
        end
        G = [G;Gij];
        H = [H;Hij];
        prob.g.complementarity(ii,jj) = {Gij.*Hij - prob.p.sigma(1), -inf, 0};
        x_prev = prob.w.x(ii,jj,n_s);
    end
end


%% Mean step equilibration Heuristic
for ii=1:N_stages
    for jj=1:N_fe
        prob.f = prob.f + prob.p.gamma_h(1)*(h0-prob.w.h(ii,jj))^2;
    end
end

%% Initialize x0
x0 = [-0.5;0.5];

%% homotopy solver
default_tol = 1e-12;
opts_casadi_nlp = struct();
opts_casadi_nlp.ipopt.print_level = 0;
opts_casadi_nlp.print_time = 0;
opts_casadi_nlp.ipopt.sb = 'yes';
opts_casadi_nlp.verbose = false;
opts_casadi_nlp.ipopt.max_iter = 500;
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
opts_casadi_nlp.ipopt.linear_solver = 'mumps';

prob.create_solver(opts_casadi_nlp);

x_res = x0;
x_curr = x0;
for step=1:N_sim
    prob.w.x(0,0,0).init = x_curr;
    prob.w.x(0,0,0).lb = x_curr;
    prob.w.x(0,0,0).ub = x_curr;
    sigma_k = 1;
    while sigma_k >= 1e-9
        prob.p.sigma(1).init = sigma_k;
        prob.solve();
        sigma_k = 0.1*sigma_k;
        prob.w.init = prob.w.res;
    end
    x_curr = prob.w.x(1,N_fe,n_s).res;
    x_res = [x_res,prob.w.x(1,1,n_s).res,prob.w.x(1,2,n_s).res,x_curr];
end

hold on
plot(x_res(1,:), x_res(2,:))
x1 = -1.7:0.001:0.7;
x2 = -(x1+0.5).^2+1;
plot(x1,x2, '--')
yline(-0.1,'--')
