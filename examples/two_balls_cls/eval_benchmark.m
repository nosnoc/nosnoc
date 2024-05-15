% eval ball benchmark
clear all;
benchmark_globals;

ref_sol_filename = "two_balls_ref_sol.mat";

%% create reference solution
[t_grid_ref, x_traj_ref, n_bounces_ref] = two_balls_spring_matlab(T_sim, x0, e, 1e-16, 0);
save(ref_sol_filename, "t_grid_ref", "x_traj_ref", "n_bounces_ref");

%% load reference solution
% load(ref_sol_filename)


% for irk_scheme = [RKSchemes.RADAU_IIA]
for irk_scheme = [RKSchemes.GAUSS_LEGENDRE]
    %%
    x_ref = x_traj_ref(end, :)';
    
    errors = cell(length(NS_VALUES), 1);
    h_values = cell(length(NS_VALUES), 1);
    labels = cell(length(NS_VALUES), 1);
    %%
    i = 1;
    for n_s = NS_VALUES 
        errors{i} = [];
        h_values{i} = [];
        
        for with_guess = [0]
        for N_sim = NSIM_VALUES
            for N_FE = NFE_VALUES
                results_filename = get_results_filename(n_s, N_sim, N_FE, irk_scheme, with_guess);
                try
                    load(results_filename);
                catch
                    disp(strcat('result not available: ', results_filename))
                    continue
                end
                if all(stats.converged)
                    x_sim_end = results.x(:, end);
                    errors{i} = [errors{i}, max(abs(x_sim_end - x_ref))];
                    h_values{i} = [h_values{i}, T_sim/(N_sim*N_FE)];
                    % h_values{i} = [h_values{i}, max(results.h)];
                    disp(strcat(results_filename, ' converged with error ', num2str(errors{i}(end), '%e')))
                else
                    disp(strcat(results_filename, ' failed.'))
                end
            end
        end
        end
        % order h
        [h_values{i}, permute] = sort(h_values{i});
        errors{i} = errors{i}(permute);
        labels{i} = get_label(settings);
        i = i+1;
    end
    
    %% order plot
    list_of_markers = {'-o','-d','-s','-x','-v','-^','-<','->','-*'};
    
    figure;
    for ii = 1:length(NS_VALUES)
        loglog(h_values{ii}, errors{ii}, list_of_markers{ii}, 'linewidth',1.5);
        hold on
    end
    

    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    xlabel('time step $\bar{h}$','interpreter','latex');
    ylabel('error $E(T)$','interpreter','latex');
    grid on
    all_h_vals = horzcat(h_values{:});
    xlim([min(all_h_vals) max(all_h_vals)])
    ylim([1e-10, 1e2])
    
    % % triangle
    % r = 7;
    % x1 = 5e-3;
    % x2 = 1e-2;
    % x3 = 1e-2;
    % y1 = 1e-7;
    % y2 = y1;
    % slope = log10(y1)+r/2;
    % y3 = 10.^(slope);
    % loglog([x1 x2], [y1 y2],'k','LineWidth',1.2)
    % loglog([x2 x3], [y2 y3],'k','LineWidth',1.2)
    % loglog([x1 x2], [y1 y3],'k','LineWidth',1.2)
    
    legend(labels, 'interpreter','latex','Location','northwest');
    
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','latex');
    filename = strcat('order_plot_', string(irk_scheme), '.pdf');
    export_fig(filename)

    fail_percentage = 100*(1 - (length(all_h_vals) / (length(NSIM_VALUES) * length(NS_VALUES))));
    disp(strcat('excluded ', num2str(fail_percentage, '%.3f'), '% of settings.'))

end


function label = get_label(settings)
    if settings.irk_scheme == "GAUSS_LEGENDRE"
        label = sprintf('Gauss-Legendre: %d', 2*settings.n_s);
    elseif settings.irk_scheme == "RADAU_IIA"
        label = sprintf('Radau IIA: %d', 2*settings.n_s-1);
        if settings.n_s == 1
            label = sprintf('implicit Euler');
        end
    else
        keyboard
    end
end
