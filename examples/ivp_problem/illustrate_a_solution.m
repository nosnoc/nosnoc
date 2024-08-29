clear all
clc
close all
import casadi.*
%% Init model and settings
model = nosnoc.model.Pss();
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
%% settings
settings.n_s = 2;
settings.N_finite_elements = 2;

% Generate Model
problem_options.N_sim = 30;
problem_options.T_sim = 0.2;
model.x0 = -0.1;
% Variable defintion
x = SX.sym('x');
model.x = x;
model.c = x;
model.S = [1;-1];
% modes of the ODE
f_11 = 1;
f_12 = 3;
model.F = [f_11 f_12];
% objective
model.f_q = x^2;
model.f_q_T = (x-5/3)^2;

integrator = nosnoc.Integrator(model, problem_options, solver_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();

%%
lambda_res = integrator.get("lambda");
mu_res = integrator.get("mu");
lambda_1 = lambda_res(1,:);
lambda_2 = lambda_res(2,:);

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

