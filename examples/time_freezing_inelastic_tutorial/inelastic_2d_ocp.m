clear all;
close all;
clc;
import casadi.*
%%
linewidth = 2.5;
[settings] = default_settings_nosnoc();  
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 1;
settings.N_homotopy = 6;
settings.opts_ipopt.ipopt.max_iter = 5e2;
settings.print_level = 3;
settings.time_freezing = 1;
settings.s_sot_max = 10;
settings.s_sot_min = 0.99;
settings.homotopy_update_rule = 'superlinear';
settings.nonsmooth_switching_fun = 0;
settings.pss_lift_step_functions = 0; 
% settings.opts_ipopt.ipopt.linear_solver = 'ma57';

%%
g = 9.81;
u_max = 10;
N_stg = 20;  
N_FE = 3; 
% Symbolic variables and bounds
q = SX.sym('q',2);
v = SX.sym('v',2); 
u = SX.sym('u');
x = [q;v];
model.T = 2;
model.N_stages = N_stg;
model.N_finite_elements  = N_FE;
model.x = x;
model.u = u;
model.e = 0;
model.mu = 0.6;
model.a_n = g;
model.x0 = [0;1.5;0;0]; 
model.f_v = [0+u;-g];
model.f_c = q(2);
model.n_dim_contact = 2;
model.g_terminal = [x-[3;0;0;0]];
model.J_tangent = [1;0];
model.lbu = -u_max;
model.ubu = u_max;
model.f_q = u'*u;
%% Call nosnoc solver
[results,stats,model,settings] = nosnoc_solver(model,settings);
% [results] = polishing_homotopy_solution(model,settings,results,stats.sigma_k);
%%
qx = results.x_opt(1,:);
qy = results.x_opt(2,:);
vx = results.x_opt(3,:);
vy = results.x_opt(4,:);
t_opt = results.x_opt(5,:);
u_opt = results.u_opt(1,:);
s_opt = results.w_opt(model.ind_sot);

figure
subplot(131)
plot(qx,qy,'LineWidth',linewidth);
xlabel('$q_1$','Interpreter','latex');
ylabel('$q_2$','Interpreter','latex');
grid on
subplot(132)
plot(t_opt,vy,'LineWidth',linewidth);
hold on
plot(t_opt,vx,'LineWidth',linewidth);
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$v$','Interpreter','latex');
legend({'$v_1$','$v_2$'},'Interpreter','latex');
subplot(133)
stairs(t_opt(1:N_FE:end),[u_opt,[nan]]','LineWidth',linewidth);
hold on
stairs(t_opt(1:N_FE:end),[s_opt',nan],'LineWidth',linewidth);
legend({'$u_1(t)$','$s(t)$'},'Interpreter','latex','location','best');
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
%% 
t_grid = results.t_grid;
f_x_fun = model.f_x_fun;
x_opt = results.x_opt;
z_opt = results.z_opt;
u_opt_extended = [];
s_opt_extended = [];
for ii = 1:N_stg
   for  jj = 1:N_FE
       u_opt_extended = [u_opt_extended,u_opt(:,ii)];
       s_opt_extended  = [s_opt_extended,s_opt(ii)];
    end
end

%%
figure
subplot(221)
plot(qx,qy,'LineWidth',linewidth);
xlabel('$q_1$','Interpreter','latex');
ylabel('$q_2$','Interpreter','latex');
grid on
subplot(223)
plot(t_opt,vy,'LineWidth',linewidth);
hold on
plot(t_opt,vx,'LineWidth',linewidth);
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$v$','Interpreter','latex');
legend({'$v_1$','$v_2$'},'Interpreter','latex');
subplot(222)
stairs(t_opt(1:N_FE:end),[u_opt,[nan]]','LineWidth',linewidth);
hold on
stairs(t_opt(1:N_FE:end),[s_opt',nan],'LineWidth',linewidth);
legend({'$u_1(t)$','$s(t)$'},'Interpreter','latex','location','best');
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
subplot(224)
if 0
    stairs(t_grid,[dtdtau],'LineWidth',linewidth);
    grid on
    hold on
    xlabel('$t$','Interpreter','latex');
    ylabel('$\frac{\mathrm{d}t}{\mathrm{d}\tau}$','Interpreter','latex');
else
    plot(t_grid,t_opt,'LineWidth',linewidth);
    grid on
    hold on
    xlabel('$t$','Interpreter','latex');
    ylabel('$t(\tau)$','Interpreter','latex');
end

%     set(gcf,'Units','inches');
%     screenposition = get(gcf,'Position');
%     set(gcf,'PaperPosition',[0 0 screenposition(3:4)*2],'PaperSize',[screenposition(3:4)*2]);
%     eval(['print -dpdf -painters ' ['ocp_example'] ])


