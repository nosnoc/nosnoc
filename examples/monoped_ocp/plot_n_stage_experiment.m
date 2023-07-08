function plot_n_stage_experiment(filename)
    close all;
    output = load(filename);
    output = output.output;
    experiments = output.experiments;
    results = output.results;
    stats = output.stats;

    experiments = cell2mat(experiments);
    stages = experiments(1:2:end, 1);

    dense_results = results(1:2:end);
    sparse_results = results(2:2:end);
    dense_stats = stats(1:2:end);
    sparse_stats = stats(2:2:end);

    dense_cpu = cellfun(@(s) s.cpu_time_total, dense_stats);
    sparse_cpu = cellfun(@(s) s.cpu_time_total, sparse_stats);

    dense_iters = cellfun(@(s) sum([s.solver_stats.iter_count]), dense_stats);
    sparse_iters = cellfun(@(s) sum([s.solver_stats.iter_count]), sparse_stats);

    line_width = 4;
    marker_size = 16;
    % plot cpu_time vs N
    fig=figure;
    hold on;
    grid on;
    plot(stages, dense_cpu, 'Xr-', 'LineWidth', line_width, 'MarkerSize', marker_size)
    plot(stages, sparse_cpu, 'Xb-', 'LineWidth', line_width, 'MarkerSize', marker_size)
    xlabel('$N$', 'interpreter', 'latex', 'FontSize', 20)
    ylabel('cpu time [s]', 'interpreter', 'latex', 'FontSize', 20)
    legend(["Stewart","Step"],  'Location', 'northwest', 'FontSize', 20)
    hold off;
    saveas(fig,'cpu_time_monoped.pdf')

    fig=figure;
    hold on;
    grid on;
    plot(stages, dense_cpu./dense_iters, 'Xr-', 'LineWidth', line_width, 'MarkerSize', marker_size)
    plot(stages, sparse_cpu./sparse_iters, 'Xb-', 'LineWidth', line_width, 'MarkerSize', marker_size)
    xlabel('$N$', 'interpreter', 'latex', 'FontSize', 20)
    ylabel('cpu time per NLP iteration [s]', 'interpreter', 'latex', 'FontSize', 20)
    legend(["Stewart","Step"], 'Location', 'northwest', 'FontSize', 20)
    hold off;
    saveas(fig,'cpu_time_per_iter_monoped.pdf')

    
end
