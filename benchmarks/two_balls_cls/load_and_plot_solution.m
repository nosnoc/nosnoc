close all
for n_s = [1]
for N_sim = [722]
    rk_scheme = RKSchemes.GAUSS_LEGENDRE;
    with_guess = 0;
    N_FE = 2;
    
    results_filename = get_results_filename(n_s, N_sim, N_FE, rk_scheme, with_guess);
    load(results_filename);
    label = strcat('n_s', num2str(n_s), ' N_sim', num2str(N_sim));
    
    if all(stats.converged)
        plot_two_ball_traj(results, label);
    else
        disp(strcat(results_filename, ' failed.'))
    end
end
end