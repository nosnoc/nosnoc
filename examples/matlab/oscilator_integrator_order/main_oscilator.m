
clear all
clc
close all

import casadi.*

% digits(64)
scenario_name = 'irk_without_fesd';
scenario_name = 'irk_fesd_legendre';
use_fesd = 1;
plot_solution_trajectory = 0;
%% Benchmark settings
% discretization settings
N_stages  = 3;
N_finite_elements = 1;

d_vec = [2 3 4 5 6];
M_vec  = round(logspace(1.5,3.4,10))+8;
% % short
% M_vec  = round(logspace(1.4,2.5,10))+7;
% 
% d_vec = [2 3 4];
% M_vec = [50 100 200];

d_vec = [2 3];
M_vec = [25 35];

legend_str = {'Implicit Euler','IRK Radau 3','IRK Radau 5','IRK Radau 7','IRK Radau 9','IRK Radau 11','IRK Radau 13'};
legend_str = {'IRK GL-2','IRK GL-3','IRK GL-4','IRK GL-5','IRK GL-6','IRK GL-7','IRK GL-8'};

legend_str = [legend_str(d_vec)];

%% settings
% collocation settings
settings = default_settings_fesd();
settings.collocation_scheme = 'radau';     % Collocation scheme: radau or legendre
% settings.collocation_scheme = 'legendre';     % Collocation scheme: radau or legendre
% settings.collocation_scheme = 'lobbato';     % Collocation scheme: radau or legendre

settings.initial_theta = 0.5;
settings.initial_lambda = 0.5;
% MPCC settings
settings.solver_name = 'solver_fesd';           % name of solver function.
settings.mpcc_mode = 3;                    % 1 - extact, 2 - smooth  ,3 -relax , 4 - penalty, 5 - l1 penalty
settings.s_elastic_max = 1e1;              % upper bound for elastic variables
% Penalty/Relaxation paraemetr
comp_tol = 1e-16;
ipopt_tol = 1e-16;
settings.comp_tol = comp_tol;
settings.sigma_0 = 1e0;                     % starting smouothing parameter
settings.sigma_N = 1e-15;                     % starting smouothing parameter
settings.N_homotopy = 30 ;% number of steps
settings.kappa = 0.05;                      % decrease rate

%^ IPOPT Settings
opts_ipopt.verbose = false;
opts_ipopt.ipopt.max_iter = 800;
opts_ipopt.ipopt.mu_strategy = 'adaptive';
opts_ipopt.ipopt.mu_oracle = 'quality-function';
opts_ipopt.ipopt.tol = ipopt_tol;
opts_ipopt.ipopt.print_level = 0;
opts_ipopt.ipopt.honor_original_bounds = 'yes';
opts_ipopt.ipopt.bound_relax_factor = 1e-16;
settings.opts_ipopt = opts_ipopt;

% finite elements with switch detection
settings.use_fesd = use_fesd;       % turn on moving finite elements algortihm
settings.fesd_complementartiy_mode = 3;       % turn on moving finite elements algortihm
settings.gamma_h = 1;                    % how much can h addapt

%% Time settings
omega = 2*pi;
% analytic solution
x_star = [exp(1);0];
s_star = [exp(2)  0; exp(2)*2*omega exp(2)];
t1_star = 1; % optimal siwtch points
T = 2;                            % time budget of transformed pseudo time
T_sim = T;
model.N_stages = N_stages;
model.N_finite_elements = N_finite_elements;

%% for results storing
errors_all_experiments = [];
errors_switch_detection_1_all_experiments = [];
errors_switch_detection_2_all_experiments = [];
complementarity_all_experiments = [];
nominal_h_all_experiments = [];
M_true_all_experiment  = [];

%% Run experiment
h_opt_full = [];

for i = 1:length(d_vec)
    d = d_vec(i);
    n_col = N_stages*(d+1); % number of collocation points per 2 finite elements
    settings.d = d; % update collocation order

    % store data for fixed d and variable M/h
    errors_current_experiment = [];
    complementarity_current_experiment = [];
    nominal_h_current_experiment = [];
    M_true_current_experiment  = [];

    for  j = 1:length(M_vec)
        M = M_vec(j);
        N_sim = round(M/n_col); % total number of simulation intevals
        h = T_sim/(N_stages*N_finite_elements*N_sim);    % nominal step lenght of a single finite elemnt;
        M_true_current = N_sim*n_col;

        % update step size
        model.T = N_stages*N_finite_elements*h;
        model.h = h;
        model.T_sim = T_sim;
        fprintf('Scheme with d = %d, current collocattion points %d , run: %d of %d \n',d,M_true_current,j,length(M_vec))

        % generate new model with updated settings;
        model = oscilator(model);
        [model,settings] = model_reformulation_fesd(model,settings);
        [results,stats] = integrator_fesd(model,settings);

         % numerical error
        x_fesd = results.x_res(:,end);
        error_x = norm(x_fesd-x_star,"inf");

        % complementarity
        max_complementarity_exp = max(stats.complementarity_stats);
        % save date current experiemnt
        errors_current_experiment = [errors_current_experiment,error_x];
        fprintf('Error with (h = %2.5f, M = %d, d = %d ) is %5.2e : \n',h,M,d,error_x);
        fprintf('Complementarity residual %5.2e : \n',max_complementarity_exp);

        complementarity_current_experiment = [complementarity_current_experiment,max_complementarity_exp];
        nominal_h_current_experiment = [nominal_h_current_experiment,h];
        M_true_current_experiment = [M_true_current_experiment,M_true_current];
    end
    errors_all_experiments = [errors_all_experiments;errors_current_experiment];
    nominal_h_all_experiments = [nominal_h_all_experiments;nominal_h_current_experiment];
    M_true_all_experiment = [M_true_all_experiment;M_true_current_experiment];
    complementarity_all_experiments = [complementarity_all_experiments;complementarity_current_experiment];
end

%% Error plots
% numerical error
figure
    for ii = 1:length(d_vec)
        loglog(M_true_all_experiment(ii,:),errors_all_experiments(ii,:),'-o','linewidth',1.5);
        hold on
    end
xlabel('$M$','interpreter','latex');
ylabel('$E(2)$','interpreter','latex');
grid on
ylim([1e-14 100])
legend(legend_str,'interpreter','latex');
saveas(gcf,[scenario_name '_error_M'])


%% some stats
if length(M_vec) == 1
    figure
    subplot(211)
    stairs(stats.homotopy_iteration_stats)
    xlabel('$N$','interpreter','latex');
    ylabel('homotopy iterations','interpreter','latex');
    grid on
    subplot(212)
    semilogy(stats.complementarity_stats+1e-20,'k.-')
    xlabel('$N$','interpreter','latex');
    ylabel('comp residual','interpreter','latex');
    grid on
end

%%
figure
loglog(M_vec,complementarity_all_experiments+1e-18,'-o','linewidth',1.5);
hold on
xlabel('$M$','interpreter','latex');
ylabel('comp residual','interpreter','latex');
grid on
legend(legend_str,'interpreter','latex');
saveas(gcf,[scenario_name '_comp_residual_M'])


%%
figure
for ii = 1:length(d_vec)
loglog(nominal_h_all_experiments(ii,:),complementarity_all_experiments(ii,:)+1e-18,'-o','linewidth',1.5);
hold on
end
xlabel('$h$','interpreter','latex');
ylabel('comp residual','interpreter','latex');
grid on
legend(legend_str,'interpreter','latex');
saveas(gcf,[scenario_name '_comp_residual_h'])

%% as function of step size
figure
for ii = 1:length(d_vec)
    loglog(nominal_h_all_experiments(ii,:),errors_all_experiments(ii,:),'-o','linewidth',1.5);
    hold on
    xlabel('$h$','interpreter','latex');
    ylabel('$E(2)$','interpreter','latex');
    grid on
end
ylim([1e-10 100])
legend(legend_str,'interpreter','latex');
saveas(gcf,[scenario_name '_error_h'])


%% plot_solution_trajectory
diff_states = results.x_res;
h_opt_full = results.h_vec;
if plot_solution_trajectory
    x1_opt = diff_states(1,:);
    x2_opt = diff_states(2,:);
    tgrid = cumsum([0;h_opt_full]);
    figure
    subplot(121)
        plot(tgrid,x1_opt,'linewidt',1.0);
        grid on
        hold on
        plot(tgrid,x2_opt,'linewidt',1.0);
        xline(1)
        xlabel('$t$','interpreter','latex');
        ylabel('$x(t)$','interpreter','latex');
        legend({'$x_1(t)$','$x_2(t)$'},'interpreter','latex');
%         xlim([1-1e-5 1+1e-5])
%         ylim([0-1e-5 0+1e-5])
    subplot(122)
    plot(x1_opt,x2_opt,'linewidt',1.8);
    hold on
    fimplicit(@(x,y) (x-0).^2+(y-0).^2-1^2, [-3 3],'r','linewidth',1.5)
    xlabel('$x_1$','interpreter','latex');
    ylabel('$x_2$','interpreter','latex');
    axis equal
    grid on
    [X,Y] = meshgrid(-3:0.35:3,-3:0.35:3);
    [U,V] = oscilator_eval(X,Y);
    quiver(X,Y,U,V,'Color',0.65*ones(3,1),'LineWidth',1.2,'AutoScaleFactor',3);
    % [X,Y] = meshgrid(-1:0.1:1,-1:0.1:1);
    % [U,V] = oscilator_eval(X,Y);
    % quiver(X,Y,U,V,'Color',0.6*ones(3,1),'LineWidth',1.2,'AutoScaleFactor',2);
    xlim([-exp(1) exp(1)]*1.01)
    ylim([-exp(1) exp(1)]*0.9)

end

%% Practial slopes
if length(M_vec) ==1
    fprintf('Error: %5.2e, h = %4.4f \n',errors_all_experiments,h)
end

%%
results.M_vec = M_vec;
results.d_vec = d_vec;
results.M_true_all_experiment =M_true_all_experiment;
results.errors_switch_detection_1_all_experiments =errors_switch_detection_1_all_experiments;
results.nominal_h_all_experiments =nominal_h_all_experiments;
results.errors_all_experiments =errors_all_experiments;
results.complementarity_all_experiments = complementarity_all_experiments;
save([scenario_name '.mat'],'results')