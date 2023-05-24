clear all;
close all;
clc;
import casadi.*
%%
linewidth = 2.5;
[settings] = NosnocOptions();  
settings.irk_scheme = IRKSchemes.GAUSS_LEGENDRE;
settings.n_s = 2;
settings.N_homotopy = 10;
settings.print_level = 3;
settings.opts_casadi_nlp.ipopt.max_iter = 1e3;
settings.dcs_mode = 'CLS';
% settings.equidistant_control_grid = 1;
settings.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
settings.cross_comp_mode = 1;
settings.friction_model = 'Conic';
settings.sigma_0 = 1e1;
settings.homotopy_update_slope = 0.2;
settings.homotopy_update_rule = 'superlinear';
settings.mpcc_mode = "Scholtes_ineq";
settings.psi_fun_type = 'STEFFENSON_ULBRICH';
%%
g = 9.81;
u_max = 10;
N_stages = 20;  
N_finite_elements = 2; 
% Symbolic variables and bounds
q = SX.sym('q', 2);
v = SX.sym('v', 2); 
u = SX.sym('u');
x = [q;v];
model.T = 2;
settings.N_stages = N_stages;
settings.N_finite_elements  = N_finite_elements;
model.x = x;
model.u = u;
model.e = 0;
model.mu = 0.7;
model.x0 = [0;1.0;0;0]; 
model.f_v = [u;-g];
model.f_c = q(2);
model.g_terminal = [x-[1;0;0;0]];
% model.f_terminal = 1e4*[x-[2;0;0;0]]'*[x-[2;0;0;0]];
model.J_tangent = [1 ;0 ];
model.D_tangent = [1 -1 ;0 0];
model.lbu = -u_max;
model.ubu = u_max;
model.f_q = u'*u;
%% Call nosnoc solver
solver = NosnocSolver(model, settings);
[results,stats] = solver.solve();

qx = results.x(1,:);
qy = results.x(2,:);
vx = results.x(3,:);
vy = results.x(4,:);
u_opt = results.u(1,:);

%%
figure
subplot(131)
plot(qx,qy,'LineWidth',linewidth);
xlabel('$q_1$','Interpreter','latex');
ylabel('$q_2$','Interpreter','latex');
grid on
subplot(132)
plot(results.t_grid,vy,'LineWidth',linewidth);
hold on
plot(results.t_grid,vx,'LineWidth',linewidth);
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$v(t)$','Interpreter','latex');
legend({'$v_1$','$v_2$'},'Interpreter','latex');
subplot(133)
stairs(results.t_grid_u,[u_opt,[nan]]','LineWidth',linewidth);
hold on
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$u(t)$','Interpreter','latex');
%     set(gcf,'Units','inches');
%     screenposition = get(gcf,'Position');
%     set(gcf,'PaperPosition',[0 0 screenposition(3:4)*2],'PaperSize',[screenposition(3:4)*2]);
%     eval(['print -dpdf -painters ' ['ocp_example'] ])


