function results = extract_results_from_solver(model,settings,results)
import casadi.*
unfold_struct(settings,'caller')
unfold_struct(model,'caller')
% Store differential states
w_opt = full(results.x);
diff_states = w_opt(ind_x);
algebraic_states = w_opt(ind_z);

u_opt = w_opt(ind_u);
u_opt = reshape(u_opt,n_u,N_stages);

if time_optimal_problem
    T_opt = w_opt(ind_t_final);
else
    T_opt = [];
end
if use_fesd
    h_opt = w_opt(ind_h);
else
    h_opt = [];
    if time_optimal_problem && ~use_speed_of_time_variables
        T = T_opt;
    end
    for ii = 1:N_stages
        h_opt = [h_opt;T/(N_stages*N_finite_elements(ii))*ones(N_finite_elements(ii),1)];
    end
end

x_opt_extended = w_opt(ind_x);
x_opt_extended  = reshape(x_opt_extended,n_x,length(x_opt_extended)/n_x);
x_opt  = x_opt_extended(:,1:n_s:end);


switch pss_mode
    case 'Stewart'
        z_opt_extended = reshape(algebraic_states,n_z,length(algebraic_states)/n_z);
        z_opt  = z_opt_extended(:,1:n_s:end);
        theta_opt_extended = [z_opt_extended(1:n_theta,:)];
        lambda_opt_extended = [z_opt_extended(n_theta+1:2*n_theta,:)];
        mu_opt_extended = [z_opt_extended(end-n_sys+1:end,:)];
        %
        theta_opt= theta_opt_extended(:,1:n_s+1:end);
        lambda_opt= lambda_opt_extended(:,1:n_s+1:end);
        mu_opt= mu_opt_extended(:,1:n_s+1:end);
    case 'Step'
        z_opt_extended = reshape(algebraic_states,n_z,length(algebraic_states)/n_z);
        z_opt  = z_opt_extended(:,1:n_s:end);
        alpha_opt_extended = [z_opt_extended(1:n_alpha,:)];
        lambda_0_opt_extended = [z_opt_extended(n_alpha+1:2*n_alpha,:)];
        lambda_1_opt_extended = [z_opt_extended(2*n_alpha+1:3*n_alpha,:)];
        %
        alpha_opt= alpha_opt_extended(:,1:n_s+1:end);
        lambda_0_opt= lambda_0_opt_extended(:,1:n_s+1:end);
        lambda_1_opt= lambda_1_opt_extended(:,1:n_s+1:end);
end
t_grid = cumsum([0;h_opt]);

%% Adapt the grid in case of time optimal problems
if time_optimal_problem
    if use_speed_of_time_variables
        s_sot = w_opt(ind_sot);
        if ~local_speed_of_time_variable
            s_sot = s_sot*ones(N_stages,1);
        end
        h_rescaled = [];
        ind_prev = 1;
        for ii = 1:N_stages
            h_rescaled = [h_rescaled;h_opt(ind_prev:N_finite_elements(ii)+ind_prev-1).*s_sot(ii)];
            ind_prev = ind_prev+N_finite_elements(ii);
        end
        t_grid = cumsum([0;h_rescaled]);
    else
        t_grid = cumsum([0;h_opt]);
    end
end

ind_t_grid_u = cumsum([1; N_finite_elements]);

results.x_opt = x_opt;
results.x_opt_extended = x_opt_extended;

results.z_opt = z_opt;
results.z_opt_extended = z_opt_extended;

results.t_grid = t_grid;
results.t_grid_u = t_grid(ind_t_grid_u);

switch pss_mode
    case 'Stewart'
        results.theta_opt = theta_opt;
        results.lambda_opt = lambda_opt;
        results.mu_opt = mu_opt;

        results.theta_opt_extended = theta_opt_extended;
        results.lambda_opt_extended = lambda_opt_extended;
        results.mu_opt_extended = mu_opt_extended;
    case 'Step'
        results.alpha_opt = alpha_opt;
        results.lambda_0_opt= lambda_0_opt;
        results.lambda_1_opt = lambda_1_opt;

        results.alpha_opt_extended = alpha_opt_extended;
        results.lambda_0_opt_extended= lambda_0_opt_extended;
        results.lambda_1_opt_extended = lambda_1_opt_extended;
end
results.u_opt = u_opt;
results.f_opt = full(results.f);
results.f = [];
results.T_opt = T_opt;
results.w_opt = w_opt;
results.h_opt = h_opt;

end
