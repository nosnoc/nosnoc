function results_filename = get_results_filename(n_s, N_sim, N_FE, irk_scheme, with_guess)
    % results_filename = strcat("two_balls_ns_", num2str(n_s), '_Nsim_', num2str(N_sim),...
    %                     "_NFE_", num2str(N_FE), '_', string(irk_scheme), '_init_', num2str(with_guess));

    results_filename = strcat("ex_", num2str(n_s), '_', num2str(N_sim),...
                        "_", num2str(N_FE), '_', string(irk_scheme), '_', num2str(with_guess));
end