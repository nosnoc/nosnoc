function [varargout] = solve_simple_car_ocp(varargin)
% the first input are the settings
% the second are optinally the model (which will complement the model function later)
import casadi.*
if nargin == 0
    [settings] = default_settings_fesd();  %% Optionally call this function to have an overview of all options.
elseif nargin == 1
    settings = varargin{1};
    unfold_struct(settings,'caller');
else
    settings = varargin{1};
    model = varargin{2};
    unfold_struct(model,'caller');
    unfold_struct(settings,'caller');
end
   
%% settings completed
settings = fill_in_missing_settings(settings);
%% Generate Model
model.time_optimal_problem = settings.time_optimal_problem;
model.time_freezing = settings.time_freezing;
model = simple_car_model(model);
%% Reformulation of the PSS into a DCS
[model,settings] = model_reformulation_fesd(model,settings);
%% Formulate NLP;
% if 0
%     [solver,solver_initalization, model,settings] = create_nlp_fesd(model,settings);
% else
    [solver,solver_initalization, model,settings] = create_nlp_fesd_develop(model,settings);
% end
%% Get variables into main workspace
unfold_struct(model,'caller');
unfold_struct(settings,'caller');
unfold_struct(solver_initalization,'caller');
%% Solve NLP
try
    [sol,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);
    solver_message = sprintf('Test sucessful with complementarity residual = %2.2e.\n',stats.complementarity_stats(end));
    results.status = 1;
catch
    solver_message = 'Test faild. Homotopy loop diverged,\n';
    results.status = 0;
end

if stats.complementarity_stats(end) > 1e-3
    solver_message = sprintf('Homotopy loop terminated with large complementarity residual: %2.2e . \n',stats.complementarity_stats(end));
    results.status = 0;
end
%% Process results
w_opt = full(sol.x);
diff_states = w_opt(ind_x);
for i = 1:n_x
    eval( ['x' num2str(i) '_opt = diff_states(' num2str(i) ':n_x+n_x*d:end);']);
end

% Check the objective
if ~isempty(ind_t_final)
    T_opt = w_opt(ind_t_final);
    clocl_state_diff = x5_opt(end) - T_opt;
end

if ~isempty(ind_t_final) && ~use_speed_of_time_variables
    T_opt = w_opt(ind_t_final);
end

if use_speed_of_time_variables
    T_opt = x5_opt(end);
end
if time_optimal_problem
    optimal_solution_message = sprintf('Optimal time is:  %2.3f \n',T_opt);
else
    optimal_solution_message = sprintf('Given numerical time is:  %2.3f \n',T);
    T_opt = T;
end
time_message = sprintf('Provided numerical time T = %2.2f final phyiscal time T_phy =  %2.2f \n',T,x5_opt(end));

if time_optimal_problem
    if abs(T_opt-(13+1/3))<1e-1
        objective_message = sprintf('Very good optimal solution found. \n');
    else
        objective_message = sprintf('Found a decent local minima. \n');
    end
else
    objective_message = sprintf('..... \n');
end
%% collect results
messages.solver_message = solver_message;
messages.optimal_solution_message = optimal_solution_message ;
messages.time_message = time_message;
messages.objective_message  = objective_message;

results.sol = sol;
results.T_opt = T_opt;
results.messages = messages;
results.w_opt = w_opt;
results.W = sol.W;
%% Outputs
varargout{1} = results;
varargout{2} = stats;
varargout{3} = solver_initalization;
varargout{4} = settings;
varargout{5} = model;
% [results,stats,solver_initalization,settings,model]
end

