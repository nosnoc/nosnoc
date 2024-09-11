% BSD 2-Clause License

% Copyright (c) 2022, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% This file is part of NOSNOC.

function [results] = integrator_order_experiment(settings,legend_str)
    import casadi.*

    %% Load user settings and model details
    unfold_struct(settings,'caller')
    %% Benchmark settings
    % discretization settings
    N_stages = 1;
    N_finite_elements = 2;
    list_of_markers = {'-o','-d','-s','-x','-v','-^','-<','->','-*'};
    %% Experiment Set Up
    % preproces of step-size  %to avoid exact switch detection by hitting the switch exactly with the given grid
    T_sim = pi/2;
    ts = 1; % eact switching time
    N_sim_vec = round(logspace(log10(N_start),log10(N_end),N_samples));
    N_sim_vec = round(N_sim_vec);
    % make all number odd
    % N_sim_vec(mod(N_sim_vec,2)==0) = N_sim_vec(mod(N_sim_vec,2)==0)+1;
    % % vector with nonminal step sizes of "outer ingeration steps"

    legend_str = [legend_str(n_s_vec)];

    %% settings
    problem_options = nosnoc.Options();
    solver_options = nosnoc.solver.Options();
    problem_options.rk_representation = rk_representation;
    problem_options.rk_scheme = rk_scheme;
    problem_options.use_fesd = use_fesd;
    problem_options.print_level = 0;
    problem_options.equidistant_control_grid = 0;
    problem_options.N_stages = N_stages;
    problem_options.N_finite_elements = N_finite_elements;
    
    % problem_options.mpcc_mode = 5;
    solver_options.homotopy_steering_strategy = 'DIRECT';
    complementarity_tol = 1e-13;
    solver_options.complementarity_tol = complementarity_tol;
    solver_options.N_homotopy = 40; % number of steps
    solver_options.homotopy_update_slope = 0.05; % decrease rate
    solver_options.use_previous_solution_as_initial_guess = 1; % warm start integrator
    solver_options.print_level = 0;
    %% Time settings
    % T = 2;                            
    omega = 2*pi;
    T = T_sim;
    x_star = [exp(T-1)*cos(2*pi*(T-1));-exp((T-1))*sin(2*pi*(T-1))];
    t1_star = 1; % optimal siwtch points

    %% for results storing
    errors_all_experiments = [];
    complementarity_all_experiments = [];
    nominal_h_all_experiments = [];
    M_true_all_experiment  = [];

    %% Run experiment
    h_opt_full = [];
    for i = 1:length(n_s_vec)
        n_s = n_s_vec(i);
        n_col = N_stages*(n_s+1); % number of collocation points per 2 finite elements
        problem_options.n_s = n_s; % update collocation order
                            % store data for fixed d and variable M/h
        errors_current_experiment = [];
        complementarity_current_experiment = [];
        nominal_h_current_experiment = [];
        M_true_current_experiment  = [];

        for  j = 1:length(N_sim_vec)
            N_sim = N_sim_vec(j);
            h_sim = T_sim/(N_sim*N_stages*N_finite_elements);
            M_true_current = N_sim*n_col;
            % update step size
            problem_options.T_sim = T_sim;
            problem_options.N_sim = N_sim;
            % generate new model with updated settings;
            model = oscilator(struct);
            integrator = nosnoc.Integrator(model, problem_options, solver_options);
            [t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();
            stats = integrator.stats;
            % numerical error
            x_fesd = x_res(:,end);
            error_x = norm(x_fesd-x_star,"inf");
            all_comp_stats = cellfun(@(x) x(end), {stats.complementarity_stats});
            max_complementarity_exp = max(all_comp_stats(end,:));
            fprintf([char(problem_options.rk_scheme) ' scheme with n_s = %d, total stage points: %d , run: %d of %d \n'],n_s,M_true_current,j,length(N_sim_vec))
            fprintf('Error with (h = %2.5f, M = %d, n_s = %d ) is %5.2e : \n',h_sim ,M_true_current,n_s,error_x);
            fprintf('Complementarity residual %5.2e : \n',max_complementarity_exp);
            % save date current experiemnt
            errors_current_experiment = [errors_current_experiment,error_x];
            complementarity_current_experiment = [complementarity_current_experiment,max_complementarity_exp];
            nominal_h_current_experiment = [nominal_h_current_experiment,h_sim];
            M_true_current_experiment = [M_true_current_experiment,M_true_current];
        end
        errors_all_experiments = [errors_all_experiments;errors_current_experiment];
        nominal_h_all_experiments = [nominal_h_all_experiments;nominal_h_current_experiment];
        M_true_all_experiment = [M_true_all_experiment;M_true_current_experiment];
        complementarity_all_experiments = [complementarity_all_experiments;complementarity_current_experiment];
    end

    %% Error plots
    figure
    for ii = 1:length(n_s_vec)
        loglog(M_true_all_experiment(ii,:),errors_all_experiments(ii,:),list_of_markers{ii},'linewidth',1.5);
        hold on
    end
    xlabel('$M$','interpreter','latex');
    ylabel('$E(2)$','interpreter','latex');
    grid on
    ylim([1e-14 100])
    legend(legend_str,'interpreter','latex');
    if save_results
        try
            saveas(gcf,['results/' scenario_name '_error_M'])
        catch
            cd ..
            saveas(gcf,['results/' scenario_name '_error_M'])
        end
    end

    %% complementarity residual
    figure
    loglog(M_true_all_experiment,complementarity_all_experiments+1e-18,list_of_markers{ii},'linewidth',1.5);
    hold on
    xlabel('$M$','interpreter','latex');
    ylabel('comp residual','interpreter','latex');
    grid on
    legend(legend_str,'interpreter','latex');
    if save_results
        saveas(gcf,['results/' scenario_name '_comp_residual_M'])
    end
    %% error as function of step size
    figure
    for ii = 1:length(n_s_vec)
        loglog(nominal_h_all_experiments(ii,:),complementarity_all_experiments(ii,:)+1e-18,list_of_markers{ii},'linewidth',1.5);
        hold on
    end
    xlabel('$h$','interpreter','latex');
    ylabel('comp residual','interpreter','latex');
    grid on
    legend(legend_str,'interpreter','latex');

    if save_results
        saveas(gcf,['results/' scenario_name '_comp_residual_h'])
    end
    %% as function of step size
    figure
    for ii = 1:length(n_s_vec)
        loglog(nominal_h_all_experiments(ii,:),errors_all_experiments(ii,:),list_of_markers{ii},'linewidth',1.5);
        hold on
        xlabel('$h$','interpreter','latex');
        ylabel('$E(2)$','interpreter','latex');
        grid on
    end
    ylim([1e-14 100])
    legend(legend_str,'interpreter','latex');
    if save_results
        saveas(gcf,['results/' scenario_name '_error_h'])
    end

    %% Practial slopes
    practical_slopes_all_experiments = [];
    for ii = 1:length(n_s_vec)
        practical_slopes_all_experiments = [practical_slopes_all_experiments;diff(log(errors_all_experiments(ii,:)))./diff(log(nominal_h_all_experiments(ii,:)))];
    end
    %%
    figure
    bar([practical_slopes_all_experiments,mean(practical_slopes_all_experiments,2)])
    xticks(1:length(n_s_vec))
    xticklabels(legend_str);
    grid on
    if save_results
        saveas(gcf,['results/' scenario_name '_practical_slopes'])
    end
    %%
    % results.M_vec = M_vec;
    results.n_s_vec = n_s_vec;
    results.M_true_all_experiment = M_true_all_experiment;
    results.nominal_h_all_experiments = nominal_h_all_experiments;
    results.legend_str = legend_str;
    results.errors_all_experiments = errors_all_experiments;
    results.n_s_vec = n_s_vec;
    results.complementarity_all_experiments = complementarity_all_experiments;
    results.practical_slopes_all_experiments = practical_slopes_all_experiments;
    if save_results
        save(['results/' scenario_name '.mat'],'results')
    end

    % close all
end

