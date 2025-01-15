clear all;
close all;
clc;
import casadi.*
%%
g = 9.81;
u_max = 10;
N_stg = 20;
N_FE = 3; 
%%
linewidth = 2.5;
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 3;
problem_options.time_freezing = 1;
problem_options.s_sot_max = 10;
problem_options.s_sot_min = 0.99;
problem_options.time_freezing_nonsmooth_switching_fun = 0;
problem_options.pss_lift_step_functions = 0;
problem_options.T = 2;
problem_options.N_stages = N_stg;
problem_options.N_finite_elements  = N_FE;
problem_options.a_n = g;

solver_options.N_homotopy = 7;
solver_options.opts_casadi_nlp.ipopt.max_iter = 1e3;
solver_options.print_level = 3;
solver_options.homotopy_update_rule = 'superlinear';
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';

%%
% Symbolic variables and bounds
q = SX.sym('q',2);
v = SX.sym('v',2); 
u = SX.sym('u');
x = [q;v];

model = nosnoc.model.Cls();
model.x = x;
model.u = u;
model.e = 0;
model.mu = 0.6;
model.x0 = [0;1.5;0;0]; 
model.f_v = [0+u;-g];
model.f_c = q(2);
model.dims.n_dim_contact = 1;
model.g_terminal = [x-[3;0;0;0]];
model.J_tangent = [1;0];
model.lbu = -u_max;
model.ubu = u_max;
model.f_q = u'*u;
%% Call nosnoc solver
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();
%%
x_res = ocp_solver.get('x');
u_opt = ocp_solver.get('u');
qx = x_res(1,:);
qy = x_res(2,:);
vx = x_res(3,:);
vy = x_res(4,:);
t_opt = x_res(5,:);
s_opt = ocp_solver.get('sot');

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
stairs(t_opt(1:N_FE:end),[s_opt,nan],'LineWidth',linewidth);
legend({'$u_1(t)$','$s(t)$'},'Interpreter','latex','location','best');
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
%% 
t_grid = ocp_solver.get_time_grid();
f_x_fun = ocp_solver.dcs.f_x_fun;
x = x_res;
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


