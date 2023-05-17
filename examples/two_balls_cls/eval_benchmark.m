% eval ball benchmark
clear all;
benchmark_globals;

ref_sol_filename = "two_balls_ref_sol.mat";

%% create reference solution
[t_grid_ref, x_traj_ref, n_bounces_ref] = two_balls_spring_matlab(T_sim, x0, e, 1e-16);
save(ref_sol_filename, "t_grid_ref", "x_traj_ref", "n_bounces_ref");

%% load reference solution
% load(ref_sol_filename)

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
    for with_guess = [0, 1]
    for N_sim = NSIM_VALUES
        for N_FE = NFE_VALUES
            results_filename = get_results_filename(n_s, N_sim, N_FE, IRK_SCHEME, with_guess);
            try
                load(results_filename);
            catch
                disp(strcat('result not available: ', results_filename))
                continue
            end
            label = strcat('n_s', num2str(n_s), ' N_sim', num2str(N_sim));
            if all(stats.converged)
                x_sim_end = results.x(:, end);
                errors{i} = [errors{i}, max(abs(x_sim_end - x_ref))];
                h_values{i} = [h_values{i}, T_sim/(N_sim*N_FE)];
                disp(strcat(results_filename, ' converged with error ', num2str(errors{i}(end), '%e')))
            else
                % keyboard
                disp(strcat(results_filename, ' failed.'))
            end
        end
    end
    end
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

set(gca,'TickLabelInterpreter','latex');
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
xlabel('$h$','interpreter','latex');
ylabel('$E(T)$','interpreter','latex');
grid on

legend(labels, 'interpreter','latex','Location','southeast');


function label = get_label(settings)
    if settings.irk_scheme == "GAUSS_LEGENDRE"
        label = sprintf('Gauss-Legendre: %d', 2*settings.n_s);
    end
end
