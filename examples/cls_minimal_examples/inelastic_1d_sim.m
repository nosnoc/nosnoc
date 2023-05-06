clear all;
clc;
import casadi.*
close all
%%
settings = NosnocOptions();
settings.irk_scheme = IRKSchemes.RADAU_IIA;
settings.n_s = 2;
settings.irk_representation = 'differential';
settings.print_level = 3;
settings.N_homotopy = 5;
settings.cross_comp_mode = 3;
settings.dcs_mode = DcsMode.CLS;
settings.multiple_solvers = 0;
settings.sigma_0 = 1;
settings.mpcc_mode = "Scholtes_ineq";
% some new verbose options for debuging
settings.print_details_if_infeasible = 0;
settings.pause_homotopy_solver_if_infeasible = 0;
settings.real_time_plot = 0;
settings.no_initial_impacts = 1;
settings.opts_ipopt.ipopt.linear_solver = 'ma97';

%%
g = 9.81;
% Symbolic variables and bounds
q = SX.sym('q',1);
v = SX.sym('v',1);
model.M = 1;
model.x = [q;v];
model.e = 0.0;
model.mu = 0;
x0 = [0.6;0];
model.x0 = x0;
model.f_v = -g;
model.f_c = q;

%% Simulation setings
N_FE = 5;
T_sim = 1.0;
N_sim = 1;

model.T_sim = T_sim;
model.N_FE = N_FE;
model.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 1;

%% Call nosnoc Integrator
[results,stats,model] = integrator_fesd(model,settings);
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
% exact impulse value
Lambda_star = v2(1)-v1(end);


%%
figure
subplot(311)
plot(t_grid,qx);
hold on
plot(tt1,q1,'k')
plot(tt2,q2,'k')
legend({'$q$ - numerical','$q$ - anlyitic'},'interpreter','latex');
xlim([0 t_grid(end)])
ylim([-0.1 x0(1)+0.1])
grid on
ylabel('$q$','interpreter','latex');
xlabel('$t$','interpreter','latex');
% axis equal
subplot(312)
plot(t_grid,vx);
hold on
plot(tt1,v1,'k')
plot(tt2,v2,'k')
plot(t_grid,vx,'b.','MarkerSize',6);
legend({'$q$ - numerical','$q$ - anlyitic'},'interpreter','latex');
ylim([-3 3])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
subplot(313)
if settings.no_initial_impacts == 1
    stem(t_grid,[nan,results.all_res.Lambda_normal_opt,nan])
else
    stem(t_grid,[results.all_res.Lambda_normal_opt,nan])
end
hold on
yline(Lambda_star,'k--')
xlim([-0.01 t_grid(end)])
ylim([-0.1 max([results.all_res.Lambda_normal_opt,Lambda_star])+1])
grid on
legend({'$\Lambda$ - numerical','$\Lambda$ - anlyitic'},'interpreter','latex');
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda$','interpreter','latex');

if N_sim == 1
    fprintf('Impulse error %2.2e \n',abs(max(results.all_res.Lambda_normal_opt))-Lambda_star)
end
