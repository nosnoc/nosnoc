clear all;
clear all;
clc;
import casadi.*
close all
%%
settings = NosnocOptions();
settings.irk_scheme = IRKSchemes.RADAU_IIA;
% settings.irk_scheme = IRKSchemes.LOBATTO_IIIA;
settings.irk_representation ="differential";
settings.n_s = 3;
settings.psi_fun_type = CFunctionType.FISCHER_BURMEISTER;
settings.print_level = 3;
settings.N_homotopy = 6;
settings.cross_comp_mode = 3;
settings.dcs_mode = DcsMode.CLS;
settings.multiple_solvers = 0;
settings.sigma_0 = 1e-6;
settings.gamma_h = 0.5;
settings.mpcc_mode = "Scholtes_ineq";

% some new verbose options for debuging
settings.print_details_if_infeasible = 0;
settings.pause_homotopy_solver_if_infeasible = 0;
settings.real_time_plot = 0;
%%
g = 9.81;
% Symbolic variables and bounds
q = SX.sym('q',1);
v = SX.sym('v',1);
model.M = 1;
model.x = [q;v];
model.e = 0;
model.mu = 0;
model.a_n = 20;
x0 = [0.5;0];
% x0 = [0;-5];
model.x0 = x0;
model.f_v = -g;
model.f_c = q;

%% Simulation setings
N_FE = 2;
T_sim = 0.5;
N_sim = 20;
model.T_sim = T_sim;
model.N_FE = N_FE;
model.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 1;

%% Call nosnoc Integrator
[results,stats,model,settings,solver] = integrator_fesd(model,settings);
%% read and plot results
unfold_struct(results,'base');
qx = x_res(1,:);
vx = x_res(2,:);


t_s = sqrt(2*x0(1)/g);
tt1 = linspace(0,t_s,100);
tt2 = linspace(t_s,T_sim,100);
v1 = x0(2)-g*tt1;
q1 = x0(1)+x0(2)*tt1-g*tt1.^2/2;
if model.e == 0
    v2 = 0*(tt2-t_s);
    q2 = 0*(tt2-t_s);
else
    v2 = -model.e*v1(end)-g*(tt2-t_s);
    q2 = q1(end)+v2(1)*(tt2-t_s)-g*(tt2-t_s).^2/2;
end
%t_opt = x_res(5,:);
t_opt = t_grid;
figure
subplot(121)
plot(t_grid,qx);
hold on
plot(tt1,q1,'k')
plot(tt2,q2,'k')
axis equal
grid on
ylabel('$q_x$','interpreter','latex');
xlabel('$t$','interpreter','latex');
% axis equal
subplot(122)
plot(t_grid,vx);
hold on
plot(t_grid,vx,'bo');
plot(tt1,v1,'k')
plot(tt2,v2,'k')
ylim([-3 3])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');

%%
% results.all_res.Lambda_normal_opt
% model.problem.print
%%
if 0 
fprintf('\t\tg\t\t lbg\t\t ubg\t\t g_sym \n')
print_casadi_matrix([results.solver_ouput.g,results.solver_initialization.lbg,results.solver_initialization.ubg,model.g])
fprintf('\t\t w\t\t lbw\t\t ubw\t\t w_sym \n')
print_casadi_matrix([results.solver_ouput.x,results.solver_initialization.lbw,results.solver_initialization.ubw,model.w])
end
