clear all;
clear all;
clc;
import casadi.*
%%
[settings] = default_settings_nosnoc();  
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 1;
settings.mpcc_mode = MpccMode.Scholtes_ineq;
settings.print_level = 2;
settings.N_homotopy = 6;
settings.use_fesd = 1;
settings.time_freezing = 1;
settings.pss_lift_step_functions= 1;
settings.impose_terminal_phyisical_time  = 1;
settings.stagewise_clock_constraint = 0;
settings.nonsmooth_switching_fun = 1;
settings.pss_lift_step_functions = 0;
%%
g = 10;
% Symbolic variables and bounds
q = SX.sym('q',3); 
v = SX.sym('v',3); 
model.e = 0;
model.mu = 0.2;
model.n_dim_contact = 3;
model.x = [q;v]; 
model.a_n = g;
model.x0 = [0;0;1;2;1;0]; 
F_ext = [1;1]*0;
% norm(F_ext)
model.f_v = [F_ext;-g];
model.f_c = q(3);
model.J_tangent = [1 0;0 1; 0 0];
%% Simulation setings
N_finite_elements = 3;
T_sim = 3;
N_sim = 20;
model.T_sim = T_sim;
model.N_FE = N_finite_elements;
model.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 0;
%% Call FESD Integrator
[results,stats,model] = integrator_fesd(model,settings);
%%
qx = results.x_res(1,:);
qy = results.x_res(2,:);
qz = results.x_res(3,:);
vx = results.x_res(4,:);
vy = results.x_res(5,:);
vz = results.x_res(6,:);
t_opt = results.x_res(7,:);
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
t_grid = results.t_grid;
lambda0 = results.lambda_0_res;
lambda1 = results.lambda_1_res;
alpha = results.alpha_res;

if settings.time_freezing_nonlinear_friction_cone
theta1 = alpha(1,:)+(1-alpha(1,:)).*(alpha(2,:));
theta2 = (1-alpha(1,:)).*(1-alpha(2,:)).*(1-alpha(3,:));
theta3  = (1-alpha(1,:)).*(1-alpha(2,:)).*(alpha(3,:));
theta = [theta1;theta2;theta3];
else
theta1 = alpha(1,:)+(1-alpha(1,:)).*(alpha(2,:));
theta2 = (1-alpha(1,:)).*(1-alpha(2,:)).*(1-alpha(3,:)).*(1-alpha(4,:));
theta3  = (1-alpha(1,:)).*(1-alpha(2,:)).*(1-alpha(3,:)).*(alpha(4,:));
theta4  = (1-alpha(1,:)).*(1-alpha(2,:)).*(alpha(3,:)).*(1-alpha(4,:));
theta5  = (1-alpha(1,:)).*(1-alpha(2,:)).*(alpha(3,:)).*(alpha(4,:));
theta = [theta1;theta2;theta3;theta4;theta5];
end

n_f = model.n_theta_step;
t_grid(1) = [];
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
% x_res = results.x_res;
% t_grid = results.t_grid;
% c_eval = [];
% for ii = 1:length(x_res)
%     c_eval = [c_eval,full(c_fun(x_res(:,ii)))];
% end
% figure
% plot(t_grid,c_eval))
% grid on

