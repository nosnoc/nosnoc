%
%    This file is part of NOSNOC.
%
%    NOSNOC -- A software for NOnSmooth Numerical Optimal Control.
%    Copyright (C) 2022 Armin Nurkanovic, Moritz Diehl (ALU Freiburg).
%
%    NOSNOC is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 3 of the License, or (at your option) any later version.
%
%    NOSNOC is distributed in the hope that it will be useful,
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
function [varargout] = nosnoc_solver(varargin)
model_unedited = varargin{1};
settings_unedited = varargin{2};
if nargin > 2
    w0 = varargin{3};
end
import casadi.*
%% Create NLP element and solve OCP with homotopy
tic
% [solver,solver_initalization,model,settings] = create_nlp_nosnoc_rv(model_unedited,settings_unedited);
    [solver,solver_initalization,model,settings] = create_nlp_nosnoc(model_unedited,settings_unedited);
solver_generating_time = toc;
if settings.print_level >=2
    fprintf('Solver generated in in %2.2f s. \n',solver_generating_time);
end
%% Check provided initial guess
if exist("w0",'var')
    if length(w0) == length(solver_initalization.w0)
        solver_initalization.w0 = w0;
    else
        fprintf('Provided user guess does not have the appropiate dimension, it should be a vector of length %d, the provided vectors has a length of %d. \n',length(solver_initalization.w0),length(w0));
        fprintf('Taking the default initial guess... \n');
    end
end

%% Solve OCP with kinematics model in time-freezing
cpu_time_presolve = 0;
if settings.time_freezing && settings.virtual_forces_kinematic_iteration
    settings_unedited.virtual_forces_in_every_mode = 1;
    [w0,cpu_time_presolve,w0_unchanged] = time_freezing_kinematics_iteration(model_unedited,settings_unedited);
end
%% Solve the discrete-time OCP
[results,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);
total_time = sum(stats.cpu_time)+cpu_time_presolve;
%% Process and store results
unfold_struct(settings,'caller');
unfold_struct(model,'caller');
results = extract_results_from_solver(model,settings,results);
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
if sum(stats.cpu_time) <60
    fprintf('Total homotopy solver time: %2.3f seconds. \n',sum(stats.cpu_time));
else
    fprintf('Total homotopy solver time: %2.3f seconds /  %2.3f minutes. \n',sum(stats.cpu_time),sum(stats.cpu_time)/60);
end
fprintf('Complementarity residual (\ell_1 norm): %2.3e.\n',complementarity_iter);

if time_optimal_problem
    T_opt = results.w_opt(model.ind_t_final);
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