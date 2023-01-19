function results = extract_results_from_solver(model,settings,results)
import casadi.*
settings_bkp = settings;    
unfold_struct(settings,'caller')
unfold_struct(model,'caller')
settings = settings_bkp;
% Store differential states
w_opt = full(results.x);
diff_states = w_opt(ind_x);
algebraic_states = w_opt(ind_z);

u_opt = w_opt([ind_u{:}]);
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
x_opt_s = [cellfun(@(x) w_opt(x), structured_ind.x, 'uni', 0)];
x_opt  = [x_opt_extended(:,1), x_opt_s{end,:}];


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

%% Get structured output (These do not contain x0)
switch pss_mode
  case 'Stewart'
    theta_opt_s = cellfun(@(theta) w_opt(theta), structured_ind.theta(1:end-(~settings.right_boundary_point_explicit),:), 'uni', 0);
    lambda_opt_s = cellfun(@(lam) w_opt(lam), structured_ind.lam, 'uni', 0);
    mu_opt_s = cellfun(@(mu) w_opt(mu), structured_ind.mu, 'uni', 0);
  case 'Step'
    alpha_opt_s = cellfun(@(alpha) w_opt(alpha), structured_ind.alpha(1:end-(~settings.right_boundary_point_explicit),:), 'uni', 0);
    lambda_n_opt_s = cellfun(@(lam) w_opt(lam), structured_ind.lambda_n, 'uni', 0);
    lambda_p_opt_s = cellfun(@(lam) w_opt(lam), structured_ind.lambda_p, 'uni', 0);
end
x_i_opt = cell(n_x, 1);
x_i_opt_flat = cell(n_x, 1);
for i = 1:n_x
    x_i_opt{i} = cellfun(@(x) x(i), x_opt_s);
    x_i_opt_flat{i} = reshape(x_i_opt{i}, prod(size(x_i_opt{i})), 1);
end

switch pss_mode
  case 'Stewart'
    theta_i_opt = cell(n_theta, 1);
    theta_i_opt_flat = cell(n_theta, 1);
    % convex multiplers
    for i = 1:n_theta
        theta_i_opt{i} = cellfun(@(t) t(i), theta_opt_s);
        theta_i_opt_flat{i} = reshape(theta_i_opt{i}, prod(size(theta_i_opt{i})), 1);
    end

    lambda_i_opt = cell(n_theta, 1);
    lambda_i_opt_flat = cell(n_theta, 1);
    % lambdas
    for i = 1:n_theta
        lambda_i_opt{i} = cellfun(@(l) l(i), lambda_opt_s);
        lambda_i_opt_flat{i} = reshape(lambda_i_opt{i}, prod(size(lambda_i_opt{i})), 1);
    end

    mu_i_opt = cell(n_sys, 1);
    mu_i_opt_flat = cell(n_sys, 1);
    % mu
    for i = 1:n_sys
        mu_i_opt{i} = cellfun(@(l) l(i), mu_opt_s);
        mu_i_opt_flat{i} = reshape(mu_i_opt{i}, prod(size(mu_i_opt{i})), 1);
    end
  case 'Step'
    alpha_i_opt = cell(n_alpha, 1);
    alpha_i_opt_flat = cell(n_alpha, 1);
    % convex multiplers
    for i = 1:n_alpha
        alpha_i_opt{i} = cellfun(@(t) t(i), alpha_opt_s);
        alpha_i_opt_flat{i} = reshape(alpha_i_opt{i}, prod(size(alpha_i_opt{i})), 1);
    end

    lambda_n_i_opt = cell(n_lambda, 1);
    lambda_n_i_opt_flat = cell(n_lambda, 1);
    % lambdas
    for i = 1:n_lambda/2
        lambda_n_i_opt{i} = cellfun(@(l) l(i), lambda_n_opt_s);
        lambda_n_i_opt_flat{i} = reshape(lambda_n_i_opt{i}, prod(size(lambda_n_i_opt{i})), 1);
    end

    lambda_p_i_opt = cell(n_lambda, 1);
    lambda_p_i_opt_flat = cell(n_lambda, 1);
    % lambdas
    for i = 1:n_lambda/2
        lambda_p_i_opt{i} = cellfun(@(l) l(i), lambda_p_opt_s);
        lambda_p_i_opt_flat{i} = reshape(lambda_p_i_opt{i}, prod(size(lambda_p_i_opt{i})), 1);
    end
end

%% Populate output
% structured output
st = struct();
st.x_i_opt = x_i_opt;
st.x_i_opt = x_i_opt_flat;
switch pss_mode
  case 'Stewart'
    st.theta_i_opt = theta_i_opt;
    st.lambda_i_opt = lambda_i_opt;
    st.mu_i_opt = mu_i_opt;
    st.theta_i_opt_flat = theta_i_opt_flat;
    st.lambda_i_opt_flat = lambda_i_opt_flat;
    st.mu_i_opt_flat = mu_i_opt_flat;
  case 'Step'
    st.alpha_i_opt = alpha_i_opt;
    st.lambda_n_i_opt = lambda_n_i_opt;
    st.lambda_p_i_opt = lambda_p_i_opt;
    st.alpha_i_opt_flat = alpha_i_opt_flat;
    st.lambda_n_i_opt_flat = lambda_n_i_opt_flat;
    st.lambda_p_i_opt_flat = lambda_p_i_opt_flat;
end
results.st = st;

% Legacy output
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
