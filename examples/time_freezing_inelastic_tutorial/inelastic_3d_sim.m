clear all;
close all;
clc;
import casadi.*
%%
g = 10;
N_finite_elements = 6;
T_sim = 3;
N_sim = 21;

%% init nosnoc
problem_options = nosnoc.Options();
integrator_options = nosnoc.integrator.Options();
solver_options = integrator_options.fesd_solver_opts;
model = nosnoc.model.Cls();
%% settings
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 2;
problem_options.use_fesd = 1;
problem_options.time_freezing = 1;
problem_options.impose_terminal_phyisical_time  = 1;
problem_options.stagewise_clock_constraint = 0;
problem_options.time_freezing_nonsmooth_switching_fun = 0;
problem_options.pss_lift_step_functions = 1;
problem_options.a_n = g;

solver_options.homotopy_steering_strategy = 'ELL_INF';
solver_options.decreasing_s_elastic_upper_bound = true;
solver_options.opts_casadi_nlp.ipopt.max_iter = 1e3;
solver_options.print_level = 3;
solver_options.N_homotopy = 6;
integrator_opts.use_previous_solution_as_initial_guess = 0;
%%
% Symbolic variables and bounds
q = SX.sym('q',3); 
v = SX.sym('v',3);
model.e = 0;
model.mu = 0.2;
model.dims.n_dim_contact = 2;
model.x = [q;v]; 
model.x0 = [0;0;1;2;1;0]; 
F_ext = [1;1]*0;
% norm(F_ext)
model.f_v = [F_ext;-g];
model.f_c = q(3);
model.J_tangent = [1 0;0 1; 0 0];
%% Simulation settings
problem_options.T_sim = T_sim;
problem_options.N_sim = N_sim;
problem_options.N_finite_elements = N_finite_elements;

%% Call FESD Integrator
integrator = nosnoc.Integrator(model, problem_options, integrator_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();
%%
qx = x_res(1,:);
qy = x_res(2,:);
qz = x_res(3,:);
vx = x_res(4,:);
vy = x_res(5,:);
vz = x_res(6,:);
t_opt = x_res(7,:);
figure
plot3(qx,qy,qz);
axis equal
xlim([-0.1 2])
ylim([-0.1 2])
zlim([-0.1 1])
grid on
xlabel('$q_x$','Interpreter','latex');
ylabel('$q_y$','Interpreter','latex');
zlabel('$q_z$','Interpreter','latex');
%%
figure
plot(t_opt,vx,'LineWidth',2);
grid on
hold on
plot(t_opt,vy);
plot(t_opt,vz);
xlabel('$t$','Interpreter','latex');
ylabel('$v$','Interpreter','latex');
legend({'$t_1^\top v$','$t_2^\top v$','$n^\top v$'},'Interpreter','latex','Location','best');


%%
lambda0 = integrator.get('lambda_n');
lambda1 = integrator.get('lambda_p');
alpha = integrator.get('alpha');

if problem_options.time_freezing_nonlinear_friction_cone
    theta1 = alpha(1,:)+(1-alpha(1,:)).*(alpha(2,:));
    theta2 = (1-alpha(1,:)).*(1-alpha(2,:)).*(1-alpha(3,:));
    theta3  = (1-alpha(1,:)).*(1-alpha(2,:)).*(alpha(3,:));
    theta = [theta1;theta2;theta3];
    n_f = 3;
else
    theta1 = alpha(1,:)+(1-alpha(1,:)).*(alpha(2,:));
    theta2 = (1-alpha(1,:)).*(1-alpha(2,:)).*(1-alpha(3,:)).*(1-alpha(4,:));
    theta3  = (1-alpha(1,:)).*(1-alpha(2,:)).*(1-alpha(3,:)).*(alpha(4,:));
    theta4  = (1-alpha(1,:)).*(1-alpha(2,:)).*(alpha(3,:)).*(1-alpha(4,:));
    theta5  = (1-alpha(1,:)).*(1-alpha(2,:)).*(alpha(3,:)).*(alpha(4,:));
    theta = [theta1;theta2;theta3;theta4;theta5];
    n_f = 5;
end

figure
for ii = 1:n_f 
    subplot(1,n_f,ii)
    plot(t_grid,theta(ii,:));
    grid on
    xlabel('$\tau$','Interpreter','latex')
    ylabel(['$\theta_' num2str(ii) '$'],'Interpreter','latex')
    ylim([-0.1 1.1])
end
%% switching functions
% c_fun = model.c_fun;
% x_res = results.x;
% t_grid = results.t_grid;
% c_eval = [];
% for ii = 1:length(x_res)
%     c_eval = [c_eval,full(c_fun(x_res(:,ii)))];
% end
% figure
% plot(t_grid,c_eval))
% grid on

