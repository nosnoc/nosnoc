clear all;
close all;

T = (11/12)*pi + sqrt(3);
N_sim = 1000;
N_fe = 2;
use_fesd = false; 
n_s = 2;
[prob, data, opts,h] = integrator(T, N_sim, N_fe, use_fesd, n_s);
orig_init = prob.w.init;
prob.p.gamma_h(1).init = 1000;
%% S I M U L A T E
x_res = data.x0;
x_res_long = data.x0;
c_ind = [0;0];
h_res = [];
lambda_res = [];
lambda_res_long = [];
x_curr = data.x0;
lambda_curr = 0;
G_res = [];
H_res = [];
eta_res = [];
for step=1:N_sim
    prob.w.init = orig_init;
    prob.w.x(0,0,data.n_s).init = x_curr;
    prob.w.x(0,0,data.n_s).lb = x_curr;
    prob.w.x(0,0,data.n_s).ub = x_curr;
    prob.w.lambda(0,0,data.n_s).init = lambda_curr;
    prob.w.lambda(0,0,data.n_s).lb = lambda_curr;
    prob.w.lambda(0,0,data.n_s).ub = lambda_curr;
    success = homotopy(prob, 1, 1e-15,0.01);
    if ~success
        disp(['Failure to converge at step=' num2str(step)])
    else
        disp(['step=' num2str(step)])
    end
    x_curr = prob.w.x(1,data.N_fe,data.n_s).res;
    lambda_curr = prob.w.lambda(1,data.N_fe,data.n_s).res;
    x_sim = prob.w.x(1,:,data.n_s).res;
    x_sim_long = prob.w.x(1,:,:).res;
    lambda_sim = prob.w.lambda(1,:,data.n_s).res;
    lambda_sim_long = prob.w.lambda(1,:,:).res;
    if opts.use_fesd
        h_sim = prob.w.h(1,:).res;
    else
        h_sim = prob.p.T(1).init/(data.N_fe)*ones(data.N_fe,1)';
    end
    h_res = [h_res,h_sim];
    x_res = [x_res,x_sim];
    x_res_long = [x_res_long,x_sim_long];
    lambda_res = [lambda_res, lambda_sim];
    lambda_res_long = [lambda_res_long, lambda_sim_long];
    G_res = [G_res, full(prob.G_fun(prob.w.res))];
    H_res = [H_res, full(prob.H_fun(prob.w.res))];
    eta_res = [eta_res, full(prob.eta_fun(prob.w.res))];
end

c_res = full(prob.c_fun(x_res));
c_res_long = full(prob.c_fun(x_res_long));
t_res = [0,cumsum(h_res)];

%% Plot
figure
hold on
plot(t_res, x_res(1,:))
plot(t_res, x_res(2,:))
hold off
if opts.use_fesd
    figure
    stairs(t_res(1:end-1),h_res)
    figure
    stairs(h_res)
end
%figure
%subplot(2,1,1)
%stairs(t_res,c_ind(1,:), '-ro')
%subplot(2,1,2);
%stairs(t_res,c_ind(2,:), '-bo')
%hold off
figure
hold on
plot(x_res(1,:), x_res(2,:))
%yline(-1, 'r--')
hold off

figure
hold on
plot(c_res_long(2:end))
plot(lambda_res_long)
hold off
figure
hold on
plot(t_res(2:end),c_res(2:end))
plot(t_res(2:end),lambda_res)
hold off

figure
plot(eta_res)
