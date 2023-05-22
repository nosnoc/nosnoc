clear all;
clc;
import casadi.*
close all

%%
benchmark_globals;

%% create model
% Symbolic variables and bounds
q = SX.sym('q',2);
v = SX.sym('v',2);


%% create reference MATLAB solution
ref_sol_filename = "two_balls_guess_sol.mat";
[t_grid_guess, x_traj_guess, n_bounces_guess, lambda_normal_guess] = two_balls_spring_matlab(1.1*T_sim, x0, e, 1e-3);
save(ref_sol_filename, "t_grid_guess", "x_traj_guess", "n_bounces_guess", "lambda_normal_guess");
% load(ref_sol_filename)
tic

%% run experiments
for kk = 1:length(IRK_SCHEME)
    for with_guess = [1 0]
        for n_s = NS_VALUES
            if kk == 3
                if n_s == 1
                    n_s = 5;
                end
            end
            for N_sim = NSIM_VALUES
                for N_FE = NFE_VALUES
                    model.M = eye(2);
                    model.x = [q;v];
                    model.e = e;
                    model.mu = 0;
                    model.x0 = x0;
                    model.f_v = [-m*g+k*(q(2)-q(1)-l);-m*g-k*(q(2)-q(1)-l)];
                    model.f_c = q(1)-R;
                    % settings
                    settings = NosnocOptions();
                    settings.irk_scheme = IRK_SCHEME{kk};
                    % settings.irk_representation = 'differential';
                    settings.n_s = n_s;
                    settings.print_level = 3;
                    settings.cross_comp_mode = 1;
                    settings.dcs_mode = DcsMode.CLS;
                    settings.multiple_solvers = 0;
                    settings.no_initial_impacts = 1;
                    settings.print_details_if_infeasible = 0;
                    settings.pause_homotopy_solver_if_infeasible = 0;
                    % settings.opts_ipopt.ipopt.linear_solver = 'ma97';
                    settings.comp_tol  = 1e-13;
                    settings.sigma_N = 1e-13;
                    settings.mpcc_mode = "elastic_ineq";
                    settings.elastic_scholtes = 1;
                    settings.sigma_0 = 1e0;
                    settings.homotopy_update_slope = 0.2;
                    settings.opts_casadi_nlp.ipopt.max_iter = 1500;
                    settings.use_previous_solution_as_initial_guess  = 1;
                    %% Simulation settings
                    model.T_sim = T_sim;
                    model.N_FE = N_FE;
                    model.N_sim = N_sim;

                    %% Call nosnoc Integrator
                    initial_guess = struct();
                    initial_guess.x_traj = x_traj_guess;
                    initial_guess.t_grid = t_grid_guess;
                    initial_guess.lambda_normal_traj = lambda_normal_guess;

                    if with_guess
                        % settings.opts_casadi_nlp.ipopt.least_square_init_duals = 'yes';
                        [results, stats, model, settings, solver] = integrator_fesd(model, settings, [], initial_guess);
                    else
                        [results, stats, model, settings, solver] = integrator_fesd(model, settings, []);
                    end

                    results_filename = get_results_filename(n_s, N_sim, N_FE, settings.irk_scheme, with_guess);
                    save(results_filename, 'results', 'stats', 'settings')

                    clear model solver
                end
            end
        end
    end
end
disp('experiment loop took')
toc


%% read and plot results
q1 = results.x(1,:);
q2 = results.x(2,:);
v1 = results.x(3,:);
v2 = results.x(4,:);
t_grid = results.t_grid;

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
ylim([-max(abs([v1,v2]))-1.0 max(abs([v1,v2]))+1])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
subplot(313)
Lambda_opt = [results.Lambda_normal];
stem(t_grid(1:N_FE:end),[Lambda_opt,nan])
hold on
xlim([-0.01 t_grid(end)])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda_{\mathrm{n}}$','interpreter','latex');

%% compare
error = norm(x_traj_guess(end,:)'-results.x(:,end));
fprintf('Numerical error %2.2e \n',error);
