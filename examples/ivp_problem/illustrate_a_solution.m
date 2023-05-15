clear all
clc
close all
import casadi.*
% settings
settings = NosnocOptions();
settings.n_s = 1;     
% Generate Model
model.T_sim = 0.2;
model.N_sim = 30;
model.N_stages = 1;
model.N_finite_elements = 2;
model.x0 = -0.1;
% Variable defintion
x = MX.sym('x');
model.x = x;
model.c = x;
model.S = [1;-1];
% modes of the ODE
f_11 = [1];
f_12 = [3];
model.F = [f_11 f_12];
% objective
model.f_q = x^2;
model.f_q_T = (x-5/3)^2;

[results,stats,model] = integrator_fesd(model,settings);

%%
t_grid = results.t_grid;
x_res = results.x;
lambda_res = results.lam;
mu_res = results.mu;
mu_res = min(x_res,-x_res);
lambda_1 = x_res - mu_res;
lambda_2 = -x_res - mu_res;

figure
subplot(121)
plot(t_grid,x_res,'k','LineWidth',1.5)
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$x(t)$','Interpreter','latex')
subplot(122)
plot(t_grid,lambda_1,'k','LineWidth',1.5)
hold on
plot(t_grid,lambda_2,'Color',0.5*ones(3,1),'LineWidth',1.5)
plot(t_grid,mu_res,'k--')
% 
% plot(t_grid,[nan,lambda_res(1,:)],'k','LineWidth',1.5)
% hold on
% plot(t_grid,[nan,lambda_res(2,:)],'Color',0.5*ones(3,1),'LineWidth',1.5)
% plot(t_grid,mu_res,'k--')
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$\lambda(t),\ \mu(t)$','Interpreter','latex')
legend({'$\lambda_1(t)$','$\lambda_2(t)$','$\mu(t)$'},'Interpreter','latex','Location','best')

