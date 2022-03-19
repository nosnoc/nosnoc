clear all
clc
close all

import casadi.*

% digits(64)
scenario_name = 'irk_without_fesd_M';
moving_finite_elements_on = 0;
plot_solution_trajectory = 1;
%% Benchmark settings
% discretization settings
N_stages  = 3;
N_finite_elements = 1;

d_vec = [1 2 3 4 5];
% d_vec = 4;
M_vec  = round(logspace(1.5,2.71,10))+13;

% 
% M_vec = [35]
% d_vec = 7;

T = 2;
T_sim = 2;
ts = 1; % switching point


%% preproces of step-size to avoid exact switch detection
% h_sim_vec = logspace(log10(0.0005),log10(0.2),10);
% N_sim_vec = M_vec./N_stages*(d)
% h_sim_vec = T./N_sim_vec;

% h_i = h_sim_vec/(N_finite_elements*N_stages);
% Ns = ts./h_i;
% ind = mod(Ns,2) ==0;
% Ns(ind) = Ns(ind)+1;
% h_sim_vec = (ts*N_finite_elements*N_stages)./Ns;

legend_str_full ={'Implicit Euler','IRK Radau 3','IRK Radau 5','IRK Radau 7','IRK Radau 9','IRK Radau 11','IRK Radau 13'};
legend_str = [legend_str_full(d_vec)];


%% settings
% complementarity and nlp tolerance
comp_tol = 1e-16;
ipopt_tol = 1e-16;

% collocation settings
settings.d = 2;                            % Degree of interpolating polynomial
settings.collocation_scheme = 'radau';     % Collocation scheme: radau or legendre
settings.initial_guess = 1;              %  1  - zeros, 2 - time interpolation , 3- sinus/cos

% MPCC settings
settings.solver_name = 'solver_mfe';           % name of solver function.
settings.mpcc_mode = 5;                    % 1 - extact, 2 - smooth  ,3 -relax , 4 - penalty, 5 - l1 penalty
settings.pointwise_or_integral = 1;        % 1 - pointwise relaxation of the comp., 0 integral relaxation of the cc.
settings.s_elastic_max = 1e2;              % upper bound for elastic variables
settings.l1_scaling = 0;
% Penalty/Relaxation paraemetr

settings.comp_tol = comp_tol;
settings.sigma_0 = 1e0;                     % starting smouothing parameter
settings.sigma_N = 1e-20;                     % starting smouothing parameter
settings.N_homotopy = 20 ;% number of steps
settings.kappa = 0.05;                      % decrease rate

% IPOPT Settings
opts_ipopt.verbose = false;
opts_ipopt.ipopt.max_iter = 800;
opts_ipopt.ipopt.mu_strategy = 'adaptive';
opts_ipopt.ipopt.mu_oracle = 'quality-function';
opts_ipopt.ipopt.tol = ipopt_tol;
opts_ipopt.ipopt.print_level = 0;
opts_ipopt.ipopt.honor_original_bounds = 'yes';
opts_ipopt.ipopt.bound_relax_factor = 0;
% opts_ipopt.ipopt.dual_inf_tol = ipopt_tol;              % Desired threshold for the dual infeasibility.
% opts_ipopt.ipopt.constr_viol_tol = ipopt_tol;           % Desired threshold for the constraint violation.
% opts_ipopt.ipopt.compl_inf_tol = ipopt_tol;           % Desired threshold for the complementarity conditions.
settings.opts_ipopt = opts_ipopt;

% finite elements with switch detection
settings.moving_finite_elements = moving_finite_elements_on;       % turn on moving finite elements algortihm
settings.gamma_h = 20;                    % how much can h addapt
settings.regularize_h =1;                 % add quadratic cost term rho_h = (h-hk)^2
settings.rho_h = 0.1;                        % regularization penalty
% settings.piecewise_equidistant_grid = 1;
%% Time settings
omega = 2*pi;
% analytic solution
x_star = [exp(1);0];
s_star = [exp(2)  0; exp(2)*2*omega exp(2)];
t1_star = 1; % optimal siwtch points
T = 2;   % time budget of transformed pseudo time
T_sim = T;

model.N_stages = N_stages;
model.N_finite_elements = N_finite_elements;
%% Vectors for results
errors_all_experiments = [];
errors_switch_detection_1_all_experiments = [];
errors_switch_detection_2_all_experiments = [];
complementarity_all_experiments = [];
nominal_h_all_experiments = [];
M_true_all_experiment  = [];
%% Run experiment

for i = 1:length(d_vec)
    d = d_vec(i);
    n_col = N_finite_elements*N_stages*(d+1); % number of stages/collocation points per integration step
    settings.d = d; % update collocation order
    % store data for fixed d and variable M/h
    errors_current_experiment = [];
    errors_switch_detection_1_current_experiment = [];
    complementarity_current_experiment = [];
    nominal_h_current_experiment = [];
    M_true_current_experiment  = [];
    for  j = 1:length(M_vec)
        N_sim = round(M_vec(j)/N_stages*(d+1));
        if mod(1/(T_sim/N_sim),2) == 0
            N_sim = N_sim +1;
        end
        h_opt_full = [];
        h_outside = T_sim/N_sim; % integrator step of FESD (outside step length)
%         N_sim = T_sim/h_outside;
        h_inside = h_outside/(N_finite_elements);    % nominal step lenght of a single finite elemnt;
        M_true_current = N_sim*n_col;

        % update step size
        model.T = h_outside;
        model.h = h_inside/N_stages;

        % generate new model with updated settings;
        model = oscilator(model);
        [model,settings] = model_reformulation_mfe(model,settings);
        [solver,solver_initalization, model,settings] = create_nlp_mfe_develop(model,settings); % generate new NLP with updated settings

        % Update Data
        unfold_struct(model,'base');
        unfold_struct(settings,'base');
        unfold_struct(solver_initalization,'base');

        sigma_k = sigma_0;
        x0 = model.x0;

        complementarity_stats = [];
        cpu_time = [];
        diff_states = [x0];
        controls = [];
        homotopy_iteration_stats = [];
        w0_base = w0;

        for jj = 1:N_sim
            % simulation for current experimetn
            tic
            solver_initalization.lbw(1:n_x) = x0;
            solver_initalization.ubw(1:n_x) = x0;
            [sol,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);
            w_opt = full(sol.x);
            complementarity_iter = full(comp_res(w_opt));
            homotopy_iteration_stats = [homotopy_iteration_stats;stats.homotopy_iterations];
            solver_initalization.w0 = w_opt;

            % read result
            x_opt = w_opt(ind_x);
            u_opt = w_opt(ind_u);
            if N_stages == 1
                x_opt = reshape(x_opt,n_x,N_finite_elements*(d+1)+1);
            else
                x_opt = reshape(x_opt,n_x,N_stages*(d+1)+1);
            end
            
            if N_stages == 1
                diff_states = [diff_states,x_opt(:,[(d+1)+1,N_finite_elements*(d+1)+1])];
            else
                diff_states = [diff_states,x_opt(:,[(d+1)+1,N_stages*(d+1)+1])];
            end

            x0 = diff_states(:,end);
            f_obj = full(sol.f);
            complementarity_stats = [complementarity_stats;complementarity_iter ];

            % stats for current step size
            if moving_finite_elements
                h_opt = w_opt(ind_h);
            else
                h_opt = h*ones(2,1);
            end
            h_opt_full = [h_opt_full ;h_opt];

        end
        fprintf('Scheme with d = %d, current collocattion points %d , run: %d of %d \n',d,M_true_current,j,length(M_vec))
        
        % numerical error
        x_mfe = diff_states(1:2,end);
        error_x = norm(x_mfe-x_star,"inf");

%         try
%             t1_mfe = tgrid(1:N_sim);
%         catch
%             t1_mfe = 1;
%         end
%         error_t1 = abs(1-t1_mfe);
        % complementarity
        max_complementarity_exp = max(complementarity_stats);

        % save date current experiemnt
        errors_current_experiment = [errors_current_experiment,error_x];
        fprintf('Error with (h = %2.5f, M = %d, d = %d ) is %5.2e : \n',h_outside,M_true_current_experiment,d,error_x);
        fprintf('Complementarity residual %5.2e : \n',max_complementarity_exp);

        complementarity_current_experiment = [complementarity_current_experiment,max_complementarity_exp];
        nominal_h_current_experiment = [nominal_h_current_experiment,h];
        M_true_current_experiment = [M_true_current_experiment,M_true_current];
        %         errors_switch_detection_1_current_experiment = [errors_switch_detection_1_current_experiment,error_t1];
    end
    % all errors (the vectors of every scheme) an n(d_vec) x n(M_vec) matrices
    % numerical error
    errors_all_experiments = [errors_all_experiments;errors_current_experiment];
    % switch detection
    %     errors_switch_detection_1_all_experiments = [errors_switch_detection_1_all_experiments;errors_switch_detection_1_current_experiment];
    nominal_h_all_experiments = [nominal_h_all_experiments;nominal_h_current_experiment];
    M_true_all_experiment = [M_true_all_experiment;M_true_current_experiment];
    complementarity_all_experiments = [complementarity_all_experiments;complementarity_current_experiment];
end


%% save results
results.M_vec = M_vec;
results.d_vec = d_vec;
results.M_true_all_experiment =M_true_all_experiment;
results.errors_switch_detection_1_all_experiments =errors_switch_detection_1_all_experiments;
results.nominal_h_all_experiments =nominal_h_all_experiments;
results.errors_all_experiments =errors_all_experiments;
results.complementarity_all_experiments = complementarity_all_experiments;
save([scenario_name '.mat'],'results')

%% Error plots
% numerical error
if 0
    figure
    loglog(M_vec,errors_all_experiments,'-o');
    xlabel('$M$','interpreter','latex');
    ylabel('$E(2)$','interpreter','latex');
    grid on
    hold on
    order_vec = d_vec*2-1;
else
    figure
    for ii = 1:length(d_vec)
        loglog(M_true_all_experiment(ii,:),errors_all_experiments(ii,:),'-o','linewidth',1.5);
        hold on
        xlabel('$M$','interpreter','latex');
        ylabel('$E(2)$','interpreter','latex');
        grid on

    end
    order_vec = d_vec*2-1;
end
%
% for ii = 1:length(order_vec)
% loglog(M_vec,+(20./(M_vec)).^order_vec(ii),'k:');
% end
ylim([1e-14 100])
legend(legend_str,'interpreter','latex');


%% some stats
if length(M_vec) == 1
    figure
    subplot(211)
    stairs(homotopy_iteration_stats)
    xlabel('$N$','interpreter','latex');
    ylabel('homotopy iterations','interpreter','latex');
    grid on
    subplot(212)
    semilogy(complementarity_stats,'k.-')
    xlabel('$N$','interpreter','latex');
    ylabel('comp residual','interpreter','latex');
    grid on
end
%%
figure
loglog(M_vec,complementarity_all_experiments,'-o','linewidth',1.5);
hold on
xlabel('$M$','interpreter','latex');
ylabel('comp residual','interpreter','latex');
grid on
legend(legend_str,'interpreter','latex');

%% as function of step size
figure
for ii = 1:length(d_vec)
    loglog(nominal_h_all_experiments(ii,:),errors_all_experiments(ii,:),'-o','linewidth',1.5);
    hold on
    xlabel('$h$','interpreter','latex');
    ylabel('$E(2)$','interpreter','latex');
    grid on
end
ylim([1e-14 100])
legend(legend_str,'interpreter','latex');

%%
% figure
% loglog(M_vec,errors_switch_detection_1_all_experiments,'-o');
% xlabel('$M$','interpreter','latex');
% ylabel('$|t_1^*-t_1|$','interpreter','latex');
% grid on
% legend(legend_str,'interpreter','latex');

%% plot_solution_trajectory
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
        xlim([1-1e-5 1+1e-5])
        ylim([0-1e-5 0+1e-5])
    %
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
% norm([x1_opt(end);x2_opt(end)]-[1;0])

%% Practial slopes
% clc
% errors_all_experiments_log = abs(log(errors_all_experiments));
% % M_vec_log = log(M_vec);
% for ii = 1:length(d_vec)
%     mean(abs(diff(errors_all_experiments_log(ii,1:3))/0.5))
% end

% if length(M_vec) ==1
% %     close all
%     fprintf('Error: %5.2e, h = %4.4f \n',errors_all_experiments,h)
% end
%%
