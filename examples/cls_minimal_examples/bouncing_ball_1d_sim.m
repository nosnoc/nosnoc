clear all;
clc;
import casadi.*
close all
%%
settings = NosnocOptions();
settings.irk_scheme = IRKSchemes.GAUSS_LEGENDRE;
settings.n_s = 2;
% settings.irk_representation = 'differential';
settings.print_level = 3;
settings.N_homotopy = 15;
settings.cross_comp_mode = 1;
settings.dcs_mode = DcsMode.CLS;
settings.multiple_solvers = 0;
settings.sigma_0 = 10;
settings.mpcc_mode = "Scholtes_ineq";
% some new verbose options for debuging
settings.print_details_if_infeasible = 0;
settings.pause_homotopy_solver_if_infeasible = 0;
settings.real_time_plot = 0;
settings.no_initial_impacts = 1;
settings.friction_model = 'Conic';
settings.conic_model_switch_handling = 'Abs';
%settings.opts_ipopt.ipopt.linear_solver = 'ma97';
settings.sigma_0 = 10;
settings.homotopy_update_slope = 0.2;

%%
g = 9.81;
% Symbolic variables and bounds
q = SX.sym('q',1);
v = SX.sym('v',1);
model.M = 1;
model.x = [q;v];
model.e = 1;
model.mu = 0;
x0 = [0.6;0];
model.x0 = x0;
model.f_v = -g;
model.f_c = q;

%% Analytic solution

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
tt = [tt1,tt2];
q_exact = [q1,q2];
v_exact = [v1,v2];

%% Simulation setings
N_FE = 6;
T_sim = 0.8;
N_sim = 1;

model.T_sim = T_sim;
model.N_FE = N_FE;
model.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 0;

%% Call nosnoc Integrator
[results,stats,model,settings,solver] = integrator_fesd(model,settings);
%% read and plot results
qx = results.x(1,:);
vx = results.x(2,:);


% exact impulse value
Lambda_star = v2(1)-v1(end);


%%
figure
subplot(311)
plot(results.t_grid,qx);
hold on
plot(tt1,q1,'k')
plot(tt2,q2,'k')
legend({'$q$ - numerical','$q$ - analytic'},'interpreter','latex');
xlim([0 results.t_grid(end)])
ylim([-0.1 x0(1)+0.1])
grid on
ylabel('$q$','interpreter','latex');
xlabel('$t$','interpreter','latex');
% axis equal
subplot(312)
plot(results.t_with_impulse, results.x_with_impulse(2,:));
hold on
plot(tt1,v1,'k')
plot(tt2,v2,'k')
plot(results.t_grid,vx,'b.','MarkerSize',6);
legend({'$q$ - numerical','$q$ - analytic'},'interpreter','latex');
ylim([-3 3])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
subplot(313)
if settings.no_initial_impacts == 1
    stem(results.t_grid,[nan,results.Lambda_normal,nan])
else
    stem(results.t_grid,[results.Lambda_normal,nan])
end
hold on
yline(Lambda_star,'k--')
xlim([-0.01 results.t_grid(end)])
ylim([-0.1 max([results.Lambda_normal,Lambda_star])+1])
grid on
legend({'$\Lambda$ - numerical','$\Lambda$ - analytic'},'interpreter','latex');
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda$','interpreter','latex');

if N_sim == 1
    fprintf('Impulse error %2.2e \n',abs(max(results.Lambda_normal)-Lambda_star))
    fprintf('position error %2.2e \n',abs(q2(end)-qx(end)))
    fprintf('velocity error %2.2e \n',abs(v2(end)-vx(end)))
end
