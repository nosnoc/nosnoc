function results = extract_results_from_solver(model, problem, settings,results)
import casadi.*
settings_bkp = settings;
unfold_struct(settings,'caller')
unfold_struct(model,'caller')
settings = settings_bkp;
% Store differential states
w_opt = full(results.nlp_results(end).x);

% populate outputs
names = {"x", "v", "z"};
switch settings.dcs_mode
  case "Stewart"
    names = [names, "theta", "lam", "mu"];
  case "Step"
    names = [names, "alpha", "lambda_n", "lambda_p"];
  case "CLS"
%     names = [names, "lambda_normal", "lambda_tangent", "y_gap", "gamma", "beta_conic", "gamma_d", "beta_d", "delta_d", "p_vt", "n_vt", "alpha_vt", "x_left_bp",...
%              "Y_gap", "Lambda_normal", "Lambda_tangent", "Gamma", "Gamma_d", "Beta_conic", "Beta_d", "Delta_d", "P_vn", "N_vn", "P_vt", "N_vt", "Alpha_vt"];
    names = [names, "lambda_normal", "lambda_tangent", "y_gap", "gamma", "beta_conic", "gamma_d", "beta_d", "delta_d", "p_vt", "n_vt", "alpha_vt", "x_left_bp",...
             "Y_gap", "Lambda_normal", "Lambda_tangent", "Gamma", "Gamma_d", "Beta_conic", "Beta_d", "Delta_d", "L_vn", "N_vt", "Alpha_vt"];
end

for name=names
    results = form_structured_output(problem, w_opt, name, results);
end

% handle x0 properly
x0 = w_opt(problem.ind_x0);
results.x_opt = [x0, results.x_opt];
results.x_opt_extended = [x0, results.x_opt_extended];




u_opt = w_opt([problem.ind_u{:}]);
u_opt = reshape(u_opt,n_u,N_stages);

if time_optimal_problem
    T_opt = w_opt(problem.ind_t_final);
else
    T_opt = [];
end
if use_fesd
    h_opt = w_opt(flatten_ind(problem.ind_h));
else
    h_opt = [];
    if time_optimal_problem && ~use_speed_of_time_variables
        T = T_opt;
    end
    for ii = 1:N_stages
        h_opt = [h_opt;T/(N_stages*N_finite_elements(ii))*ones(N_finite_elements(ii),1)];
    end
end

t_grid = cumsum([0;h_opt]);

%% Adapt the grid in case of time optimal problems
if time_optimal_problem
    if use_speed_of_time_variables
        s_sot = w_opt(flatten_ind(problem.ind_sot));
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

results.t_grid = t_grid;
results.t_grid_u = t_grid(ind_t_grid_u);


results.u_opt = u_opt;
results.f = full(results.nlp_results(end).f);
results.T_opt = T_opt;
results.w_opt = w_opt;
results.h_opt = h_opt;

end
