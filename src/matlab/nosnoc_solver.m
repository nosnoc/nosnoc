%
%    This file is part of NOSNOC.
%
%    NOS-NOC -- A software for NOnSmooth Numerical Optimal Control.
%    Copyright (C) 2022 Armin Nurkanovic, Moritz Diehl (ALU Freiburg).
%
%    NOS-NOC is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 3 of the License, or (at your option) any later version.
%
%    NOS-NOC is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with NOS-NOC; if not, write to the Free Software Foundation,
%    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%
%% Possible outputs
% [results]= nosnoc_solver(model,settings);
% [results,stats] = nosnoc_solver(model,settings);
% [results,stats,model] = nosnoc_solver(model,settings);
% [results,stats,model,settings] = nosnoc_solver(model,settings);
% [results,stats,model,settings,solver_initalization] = nosnoc_solver(model,settings);
function [varargout] = nosnoc_solver(model,settings)
import casadi.*
%% Create NLP element and solve OCP with homotopy
[solver,solver_initalization, model,settings] = create_nlp_nosnoc(model,settings);
[results,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);
total_time = sum(stats.cpu_time);
%% Process and store results
unfold_struct(settings,'caller');
unfold_struct(model,'caller');
results = extract_results_from_solver(model,settings,results);
%% heuristic polishing step for active set projection (not recomended to use yet)
if polishing_step
    [results] = polishing_homotopy_solution(model,settings,solver,stats,solver_initalization,results);
end
complementarity_iter = full(comp_res(results.w_opt));

%% Verbose
stats.total_time  = total_time;
fprintf('\n');
fprintf('-----------------------------------------------------------------------------------------------\n');
if use_fesd
    fprintf( ['OCP with the FESD ' irk_scheme ' in ' irk_representation ' mode with %d RK-stages, %d finite elements and %d control intervals.\n'],n_s,N_finite_elements(1),N_stages);
else
    fprintf( ['OCP with the Std ' irk_scheme ' in ' irk_representation ' mode with %d RK-stages, %d finite elements and %d control intervals.\n'],n_s,N_finite_elements(1),N_stages);
end
fprintf('Max homotopy iteration time: %2.3f seconds, min homotopy iteration time: %2.3f seconds.\n',max(stats.cpu_time),min(stats.cpu_time));
fprintf('Total homotopy iterations: %d.\n',stats.homotopy_iterations);
fprintf('Total homotopy solver time: %2.3f seconds. \n',sum(stats.cpu_time));
fprintf('Complementarity residual: %2.3e.\n',complementarity_iter);

if time_optimal_problem
    T_opt = results.w_opt(model.ind_t_final);
    if T_opt  < 1e-3
        T_opt = t_grid(end);
    end
    fprintf('Final time T_opt: %2.4f.\n',T_opt);
end
fprintf('-----------------------------------------------------------------------------------------------\n\n');

%% Output
varargout{1} = results;
varargout{2} = stats;
varargout{3} = model;
varargout{4} = settings;
varargout{5} = solver_initalization;
end