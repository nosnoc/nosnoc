% BSD 2-Clause License

% Copyright (c) 2022, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% This file is part of NOSNOC.

%
%
% TODO: remove varargout / varargin and make API clear!
% Examples of calling the function
% [results] = integrator_fesd(model,settings);
% [results,stats] = integrator_fesd(model,settings);
% [results,stats,model] = integrator_fesd(model,settings);

function [varargout] = integrator_fesd(model, settings, u_sim, initial_guess)

%% generate time-freezing model before turning off time-related settings
if settings.time_freezing
    [model,settings] = time_freezing_reformulation(model,settings);
end

%% Settings of integration
[model] = refine_model_integrator(model,settings);

%% Create solver functions for integrator step
solver = NosnocSolver(model, settings);
model = solver.model;
settings = solver.settings;
dims = model.dims;

unfold_struct(settings,'caller')
unfold_struct(model,'caller')

%% check does the provided u_sim has correct dims
if dims.n_u > 0
    if ~exist('u_sim','var')
        warning('no control values (u_sim) provided, using zeros')
        u_sim = zeros(n_u, N_sim);
    else
        [n_rows, n_cols] = size(u_sim);
        if n_rows ~= dims.n_u || n_cols~=N_sim
            error('Control inputs u_sim has the wrong size. Required dimension is n_u x N_sim.')
        end
    end
end

%% Initialization
results = struct();
% states
complementarity_stats  = [];
homotopy_iteration_stats = [];
time_per_iter = [];
simulation_time_pased = 0;
W = [];
sim_step_solver_results = [];
t_current = 0;
converged = [];

if exist('initial_guess', 'var')
    % remove duplicate indices
    [~, w] = unique(initial_guess.t_grid, 'stable');
    initial_guess.t_grid = initial_guess.t_grid(w);
    initial_guess.x_traj = initial_guess.x_traj(w, :);
end


%% Main simulation loop
for ii = 1:model.N_sim
    if dims.n_u > 0
        solver.set('u', {u_sim(:,ii)});
    end

    %% 
    % TODO Set up homotopy solver to take p_val explicitly
    if ii > 1 && settings.use_previous_solution_as_initial_guess
        % TODO make this possible via solver interface directly
        solver.problem.w0(n_x+1:end) = res.w(n_x+1:end);
    end

    % set all x values to xcurrent, as this is the best available guess.
    solver.set('x', {x0})
    solver.set('x_left_bp', {x0})
    if exist('initial_guess', 'var')
        t_guess = t_current + cumsum([0; model.h_k * ones(model.N_finite_elements, 1)]);
        x_guess = interp1(initial_guess.t_grid, initial_guess.x_traj, t_guess);
        % 
        x_init = cell(1, dims.N_finite_elements);
        y_gap_init = cell(1, dims.N_finite_elements);
        for j = 1:dims.N_finite_elements
            x_init{j} = x_guess(j+1, :);
            y_gap_init{j} = full(model.f_c_fun(x_init{j}));
        end
        solver.set('x', x_init);
        solver.set('x_left_bp', x_init);
        %
        if isequal(settings.dcs_mode, 'CLS')
            ind_v = dims.n_q + 1: 2*dims.n_q;
            solver.set('lambda_normal', {0});
            L_vn_init = {(x_init{1}(ind_v) + model.e * x_init{end}(ind_v))};
            solver.set('L_vn', L_vn_init);
            % 
            diff_x = x_guess(1, :) - x_guess(end, :);
            diff_v = diff_x(ind_v);
            Lambda_normal_init = max(abs(diff_v));
            % Note: this is a bit hacky..
            if Lambda_normal_init < 4
                Lambda_normal_init = 0.0;
            end
            solver.set('Lambda_normal', {Lambda_normal_init});
            disp(['init Lambda_normal', num2str(Lambda_normal_init)]);
            solver.set('y_gap', y_gap_init);
        end
    end

    %% solve
    [sol,stats] = solver.solve();
    [res, names] = extract_results_from_solver(model, solver.problem, settings, sol);
    names = [names, {"h"}];
    sim_step_solver_results = [sim_step_solver_results,res];
    time_per_iter = [time_per_iter; stats.cpu_time_total];

    % Initialize results
    if ii == 1
        for name=names
            if isfield(res, name)
                results.(name) = [];
                results.extended.(name) = [];
            end
        end
        results.x = [x0];
        results.extended.x = [x0];
        results.s_sot = [];
        results.x_with_impulse = [x0];
        results.t_with_impulse = [0];
    end

    if stats.converged == 0
        % TODO: return some infeasibility status
        warning(['integrator_fesd: did not converge in step ', num2str(ii)])
        converged = [converged, stats.converged];
    elseif print_level >=2
        fprintf('Integration step %d / %d (%2.3f s / %2.3f s) converged in %2.3f s. \n',...
            ii, model.N_sim,simulation_time_pased, model.T_sim, time_per_iter(end));
        converged = [converged, stats.converged];
    end
    simulation_time_pased = simulation_time_pased + model.T;
    

    %% update initial guess and inital value
    x0 = res.x(:,end);
    %     update clock state
    if impose_terminal_phyisical_time
        solver.problem.p0(end) = solver.problem.p0(end)+model.T;
    end
    solver.set("x0", x0);

    % TODO Set up homotopy solver to take p_val explicitly
    if use_previous_solution_as_initial_guess
        % TODO make this possible via solver interface directly
        solver.problem.w0(n_x+1:end) = res.w(n_x+1:end);
    end

    % set all x values to xcurrent, as this is the best available guess.
    solver.set('x', {x0})
    solver.set('x_left_bp', {x0})
    % try zeros?
    % solver.set('lambda_normal', 1.0)

    %% Store data
    % update results struct
    for name=names
        if name == 'x' || ~isfield(res, name)
            continue
        end
        
        results.(name) = [results.(name), res.(name)];
        if ~strcmp(name, 'h')
            results.extended.(name) = [results.extended.(name), res.extended.(name)];
        end
    end
    results.x = [results.x, res.x(:, 2:end)];
    results.extended.x = [results.extended.x, res.extended.x(:, 2:end)];

    %sot TODO: is there a better way to do this
    results.s_sot  = [results.s_sot, res.w(flatten_ind(solver.problem.ind_sot))];

    if settings.dcs_mode == DcsMode.CLS
        results.x_with_impulse = [results.x_with_impulse, res.x_with_impulse(:,2:end)];

        % TODO maybe make solver take t0 and T_final as param?
        
        results.t_with_impulse = [results.t_with_impulse, t_current + res.t_with_impulse(2:end)];
    end

    t_current = t_current + model.T;
    %stats
    complementarity_stats  = [complementarity_stats; stats.complementarity_stats(end)];
    homotopy_iteration_stats = [homotopy_iteration_stats;stats.homotopy_iterations];

    %% plot during execution
    if real_time_plot
        figure(100)
        clf
        grid on
        hold on
        xlim([0 T_sim]);
        ylabel('$x(t)$','Interpreter','latex');
        if settings.time_freezing
            plot(res.x(end,:),res.x(1:end-1,:));
            xlabel('$t$ [phyisical time]','Interpreter','latex');
            ylim([min(min(res.x(1:end-1,:)))-0.3 max(max(res.x(1:end-1,:)))+0.3])
        else
            t_temp = [0,cumsum(h_vec)'];
            plot(t_temp,res.x(:,1:end));
            xlabel('$t$','Interpreter','latex');
            ylim([min(res.x(:))-0.3 max(res.x(:))+0.3])
        end
    end
end
total_time = sum(time_per_iter);
%% Verbose
fprintf('\n');
fprintf('----------------------------------------------------------------------------------------------------------------------\n');
if use_fesd
    fprintf( ['Simulation with the FESD ' char(irk_scheme) ' with %d-RK stages completed.\n'],n_s);
else
    fprintf( ['Simulation with the standard ' char(irk_scheme) ' with %d-RK stages completed.\n'],n_s);
end
fprintf( ['RK representation: ' char(irk_representation) '.\n']);
fprintf('---------------------------------------------- Stats summary----------------------------------------------------------\n');
fprintf('N_sim\t step-size\t\tN_stg\tN_FE\t CPU Time (s)\t Max. CPU (s)/iter\tMin. CPU (s)/iter\tMax. comp.\tMin. comp.\n');
fprintf('%d\t\t\t%2.3f\t\t%d\t\t%d\t\t%2.3f\t\t\t\t%2.3f\t\t\t%2.3f\t\t\t\t%2.2e\t%2.2e\n',N_sim,h_sim,N_stages,N_finite_elements(1),total_time,max(time_per_iter),min(time_per_iter),max(complementarity_stats),min(complementarity_stats));
fprintf('----------------------------------------------------------------------------------------------------------------------\n\n');
%% Output

stats.complementarity_stats = complementarity_stats;
stats.time_per_iter = time_per_iter;
stats.homotopy_iteration_stats = homotopy_iteration_stats;
stats.converged = converged;

results.t_grid = cumsum([0,results.h])';

results.solver_ouput = sol;
results.sim_step_solver_results = sim_step_solver_results;

varargout{1} = results;
varargout{2} = stats;
varargout{3} = model;
varargout{4} = settings;
varargout{5} = solver;
end

