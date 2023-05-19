function results = extract_results_from_solver(model, problem, settings, results)
import casadi.*
% Store differential states
w_opt = full(results.nlp_results(end).x);
results.w = w_opt;

names = get_result_names_from_settings(settings);
% populate outputs
for name=names
    results = form_structured_output(problem, w_opt, name, results);
end

% handle x0 properly
x0 = w_opt(problem.ind_x0);
results.x = [x0, results.x];
results.extended.x = [x0, results.extended.x];

u = w_opt([problem.ind_u{:}]);
u = reshape(u,model.dims.n_u,model.dims.N_stages);

results.u = u;

if settings.time_optimal_problem
    T_opt = w_opt(problem.ind_t_final);
else
    T_opt = problem.model.T;
end
results.T = T_opt;

if settings.use_fesd
    h_opt = w_opt(flatten_ind(problem.ind_h))';
else
    h_opt = [];
    if settings.time_optimal_problem && ~settings.use_speed_of_time_variables
        T = T_opt;
    end
    for ii = 1:model.dims.N_stages
        h_opt = [h_opt,model.T/(model.dims.N_stages*model.dims.N_finite_elements(ii))*ones(1, model.dims.N_finite_elements(ii))];
    end
end
results.h = h_opt;

t_grid = cumsum([0,h_opt]);

%% Adapt the grid in case of time optimal problems
if settings.time_optimal_problem
    if settings.use_speed_of_time_variables
        s_sot = w_opt(flatten_ind(problem.ind_sot));
        if ~local_speed_of_time_variable
            s_sot = s_sot*ones(model.dims.N_stages,1);
        end
        h_rescaled = [];
        ind_prev = 1;
        for ii = 1:model.dims.N_stages
            h_rescaled = [h_rescaled,h_opt(ind_prev:model.dims.N_finite_elements(ii)+ind_prev-1).*s_sot(ii)];
            ind_prev = ind_prev+model.dims.N_finite_elements(ii);
        end
        t_grid = cumsum([0,h_rescaled]);
    else
        t_grid = cumsum([0,h_opt]);
    end
end
ind_t_grid_u = cumsum([1; model.dims.N_finite_elements]);

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
