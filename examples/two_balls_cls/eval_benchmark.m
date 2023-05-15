% eval ball benchmark
clear all;
benchmark_globals;


ref_sol_filename = "two_balls_ref_sol.mat";

%% create reference solution
% [t_grid_ref, x_traj_ref, n_bounces_ref] = two_balls_spring_matlab(T_sim, x0, e, 1e-12);
% save(ref_sol_filename, "t_grid_ref", "x_traj_ref", "n_bounces_ref");

%% load reference solution
load(ref_sol_filename)

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
    labels{i} = strcat('$n_s= $', num2str(n_s));
    for N_sim = NSIM_VALUES
        results_filename = strcat("two_balls_ns_", num2str(n_s), '_Nsim_', num2str(N_sim));
        load(results_filename);
        label = strcat('n_s', num2str(n_s), ' N_sim', num2str(N_sim));
        if all(stats.converged)
            disp(strcat(label, ' converged.'))
            x_sim_end = results.x_res(:, end);
            errors{i} = [errors{i}, max(abs(x_sim_end - x_ref))];
            h_values{i} = [h_values{i}, T_sim/N_sim];
        else
            disp(strcat(label, ' failed.'))
        end
    end
    i = i+1;
end


%% order plot
list_of_markers = {'-o','-d','-s','-x','-v','-^','-<','->','-*'};

figure;
for ii = 1:length(NS_VALUES)
    loglog(h_values{ii}, errors{ii}, list_of_markers{ii}, 'linewidth',1.5);
    hold on
end

set(gca,'TickLabelInterpreter','latex')
xlabel('$h$','interpreter','latex');
ylabel('$E(T)$','interpreter','latex');
grid on

legend(labels, 'interpreter','latex','Location','southwest');
