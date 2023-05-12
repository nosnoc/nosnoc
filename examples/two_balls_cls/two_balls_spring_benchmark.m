clear all;
clc;
import casadi.*
close all


NS_VALUES = [1, 2];
NSIM_VALUES = 80 * [1, 2, 4, 8, 10];
T_sim = 1;

%% create model
g = 9.81;
R = 0.2;
k = 1e4;
l = 1;
m = 1;
% Symbolic variables and bounds
q = SX.sym('q',2);
v = SX.sym('v',2);
x0 = [1;2;0;0];
e = 0.8;


%% create reference MATLAB solution
ref_sol_filename = "two_balls_ref_sol.mat";
% [t_grid_ref, x_traj_ref, n_bounces_ref] = two_balls_spring_matlab(T_sim, x0, e, 1e-11);
% save(ref_sol_filename, "t_grid_ref", "x_traj_ref", "n_bounces_ref");
load(ref_sol_filename)

%% run experiments
N_FE = 2;

for n_s = NS_VALUES
    for N_sim = NSIM_VALUES
        results_filename = strcat("two_balls_ns_", num2str(n_s), '_Nsim_', num2str(N_sim));

        %
        model.M = eye(2);
        model.x = [q;v];
        model.e = e;
        model.mu = 0;
        model.x0 = x0;
        model.f_v = [-m*g+k*(q(2)-q(1)-l);-m*g-k*(q(2)-q(1)-l)];
        model.f_c = q(1)-R;
        % settings
        settings = NosnocOptions();
        settings.irk_scheme = IRKSchemes.GAUSS_LEGENDRE;
        % settings.irk_representation = 'differential';
        settings.n_s = n_s;
        settings.print_level = 3;
        % settings.N_homotopy = 8;
        settings.cross_comp_mode = 3;
        settings.dcs_mode = DcsMode.CLS;
        settings.multiple_solvers = 0;
        settings.mpcc_mode = "Scholtes_ineq";
        settings.no_initial_impacts = 1;
        settings.sigma_0 = 1e-1;
        settings.homotopy_update_slope = 0.1;

        %% Simulation settings        
        model.T_sim = T_sim;
        model.N_FE = N_FE;
        model.N_sim = N_sim;

        %% Call nosnoc Integrator
        initial_guess = struct();
        initial_guess.x_traj = x_traj_ref;
        initial_guess.t_grid = t_grid_ref;
        settings.sigma_0 = 1e-3;
        
        [results, stats, model, settings, solver] = integrator_fesd(model, settings, [], initial_guess);

        save(results_filename, "results", "stats", "settings")

        clear model solver
    end
end


%% read and plot results
unfold_struct(results,'base');
q1 = x_res(1,:);
q2 = x_res(2,:);
v1 = x_res(3,:);
v2 = x_res(4,:);

%%
figure
subplot(311)
plot(t_grid,q1,'LineWidth',1.5);
hold on
plot(t_grid,q2,'LineWidth',1.5);
yline(R,'k--')
xlim([0 t_grid(end)])
% ylim([-1.0 max([q1,q2])+1])
grid on
ylabel('$q$','interpreter','latex');
xlabel('$t$','interpreter','latex');
% axis equal
subplot(312)
plot(t_grid,v1,'LineWidth',1.5);
hold on
plot(t_grid,v2,'LineWidth',1.5);
xlim([0 t_grid(end)])
% ylim([-max(abs([v1,v2]))-1.0 max(abs([v1,v2]))+1])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
subplot(313)
Lambda_opt = [results.all_res.Lambda_normal_opt];
stem(t_grid(1:N_FE:end),[Lambda_opt,nan])
hold on
xlim([-0.01 t_grid(end)])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda_{\mathrm{n}}$','interpreter','latex');

%% compare
error = norm(x_traj_ref(end,:)'-x_res(:,end));
fprintf('Numerical error %2.2e \n',error);



