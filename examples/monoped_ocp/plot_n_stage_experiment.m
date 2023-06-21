function plot_n_stage_experiment(filename)
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
    
    % plot cpu_time vs N
    figure;
    hold on;
    grid on;
    plot(stages, dense_cpu, 'Xr-')
    plot(stages, sparse_cpu, 'Xb-')
    xlabel('$N$', 'interpreter', 'latex')
    ylabel('cpu time [s]', 'interpreter', 'latex')
    hold off;

    figure;
    hold on;
    grid on;
    plot(stages, dense_cpu./dense_iters, 'Xr-')
    plot(stages, sparse_cpu./sparse_iters, 'Xb-')
    xlabel('$N$', 'interpreter', 'latex')
    ylabel('cpu time per NLP iteration [s]', 'interpreter', 'latex')
    hold off;
end
