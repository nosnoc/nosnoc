function [results, names] = extract_results_from_solver(model, problem, settings,results)
import casadi.*
settings_bkp = settings;
unfold_struct(settings,'caller')
unfold_struct(model,'caller')
settings = settings_bkp;
% Store differential states
w_opt = full(results.nlp_results(end).x);
results.w = w_opt;
% populate outputs
names = {"x", "v", "z"};
switch settings.dcs_mode
  case "Stewart"
    names = [names, "theta", "lam", "mu"];
  case "Step"
    names = [names, "alpha", "lambda_n", "lambda_p"];
  case "CLS"
    names = [names, "lambda_normal", "lambda_tangent", "y_gap", "gamma", "beta_conic", "gamma_d", "beta_d", "delta_d", "p_vt", "n_vt", "alpha_vt", "x_left_bp",...
             "Y_gap", "Lambda_normal", "Lambda_tangent", "Gamma", "Gamma_d", "Beta_conic", "Beta_d", "Delta_d", "L_vn", "N_vt", "Alpha_vt"];
end

for name=names
    results = form_structured_output(problem, w_opt, name, results);
end

% handle x0 properly
x0 = w_opt(problem.ind_x0);
results.x = [x0, results.x];

u = w_opt([problem.ind_u{:}]);
u = reshape(u,n_u,N_stages);

results.u = u;

if time_optimal_problem
    T_opt = w_opt(problem.ind_t_final);
else
    T_opt = problem.model.T;
end
results.T = T_opt;

if use_fesd
    h_opt = w_opt(flatten_ind(problem.ind_h))';
else
    h_opt = [];
    if time_optimal_problem && ~use_speed_of_time_variables
        T = T_opt;
    end
    for ii = 1:N_stages
        h_opt = [h_opt,T/(N_stages*N_finite_elements(ii))*ones(N_finite_elements(ii),1)];
    end
end
results.h = h_opt;

t_grid = cumsum([0,h_opt]);

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

if settings.dcs_mode == DcsMode.CLS
    x_with_impulse = [x0];
    t_with_impulse = kron(t_grid, ones(2,1));
    for ii=1:size(results.structured.x,1)
        for jj=1:size(results.structured.x,2)
            x_with_impulse = [x_with_impulse,results.structured.x_left_bp{ii,jj}];
            x_with_impulse = [x_with_impulse,results.structured.x{ii,jj}];
        end
    end
    results.x_with_impulse = x_with_impulse;
    results.t_with_impulse = t_with_impulse(1:end-1);
end

results.t_grid = t_grid;
results.t_grid_u = t_grid(ind_t_grid_u);

results.f = full(results.nlp_results(end).f);
results.g = full(results.nlp_results(end).g);
end
