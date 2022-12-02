% Copyright 2022 Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% This file is part of NOSNOC.

% The 2-Clause BSD License

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
%
%% Possible outputs
% [results]= nosnoc_solver(model,settings);
% [results,stats] = nosnoc_solver(model,settings);
% [results,stats,model] = nosnoc_solver(model,settings);
% [results,stats,model,settings] = nosnoc_solver(model,settings);
% [results,stats,model,settings,solver_initialization] = nosnoc_solver(model,settings);
function [varargout] = nosnoc_solver(varargin)
model_unedited = varargin{1};
settings_unedited = varargin{2};
if nargin > 2
    w0 = varargin{3};
end
import casadi.*
%% Create NLP element and solve OCP with homotopy
tic
% [solver,solver_initialization,model,settings] = create_nlp_nosnoc_rv(model_unedited,settings_unedited);
    [solver,solver_initialization,model,settings] = create_nlp_nosnoc(model_unedited,settings_unedited);
solver_generating_time = toc;
if settings.print_level >=2
    fprintf('Solver generated in in %2.2f s. \n',solver_generating_time);
end
%% Check provided initial guess
if exist("w0",'var')
    if length(w0) == length(solver_initialization.w0)
        solver_initialization.w0 = w0;
    else
        fprintf('Provided user guess does not have the appropiate dimension, it should be a vector of length %d, the provided vectors has a length of %d. \n',length(solver_initialization.w0),length(w0));
        fprintf('Taking the default initial guess... \n');
    end
end

%% Solve the discrete-time OCP
[results,stats,solver_initialization] = homotopy_solver(solver,model,settings,solver_initialization);
total_time = sum(stats.cpu_time);
%% Process and store results
unfold_struct(settings,'caller');
unfold_struct(model,'caller');
results = extract_results_from_solver(model,settings,results);
complementarity_iter_ell_1 = full(comp_res(results.w_opt,[model.p_val';full(model.lambda00_fun(x0))]));
switch pss_mode
    case 'Step'
    temp = [results.alpha_opt_extended.*results.lambda_0_opt_extended,(1-results.alpha_opt_extended).*results.lambda_1_opt_extended];
    complementarity_iter_ell_inf = max(temp(:));
    case 'Stewart'
        % TODO: considert cross comps as well in the inf norm
     temp = [results.theta_opt_extended.*results.lambda_opt_extended];
    complementarity_iter_ell_inf = max(temp(:));
end

%% Verbose
stats.total_time  = total_time;
fprintf('\n');
fprintf('-----------------------------------------------------------------------------------------------\n');
if use_fesd
    fprintf( ['OCP with the FESD ' irk_scheme ' in ' irk_representation ' mode with %d RK-stages, %d finite elements and %d control intervals.\n'],n_s,N_finite_elements(1),N_stages);
else
    fprintf( ['OCP with the Std ' irk_scheme ' in ' irk_representation ' mode with %d RK-stages, %d finite elements and %d control intervals.\n'],n_s,N_finite_elements(1),N_stages);
end
% fprintf('Total homotopy iterations: %d.\n',stats.homotopy_iterations);
% if sum(stats.cpu_time) <60
%     fprintf('Total homotopy solver time: %2.3f seconds. \n',sum(stats.cpu_time));
% else
%     fprintf('Total homotopy solver time: %2.3f seconds /  %2.3f minutes. \n',sum(stats.cpu_time),sum(stats.cpu_time)/60);
% end
% fprintf('Max homotopy iteration time: %2.3f seconds. \nMin homotopy iteration time: %2.3f seconds.\n',max(stats.cpu_time),min(stats.cpu_time));
% fprintf('Complementarity residual (1-norm): %2.2e.\n',complementarity_iter_ell_1);
% fprintf('Complementarity residual (inf-norm): %2.2e.\n',complementarity_iter_ell_inf);
%%
fprintf('---------------------------------------------- Stats summary----------------------------------------------------------\n');
if sum(stats.cpu_time) < 60
    fprintf('H. iters\t CPU Time (s)\t Max. CPU (s)/iter\tMin. CPU (s)/iter\t Comp. res. (ell_1) \tComp. res. (ell_inf)\n');
    fprintf('%d\t\t\t\t%2.2f\t\t\t%2.2f\t\t\t\t%2.2f\t\t\t\t%2.2e\t\t\t\t%2.2e\n',stats.homotopy_iterations,sum(stats.cpu_time),max(stats.cpu_time),min(stats.cpu_time),complementarity_iter_ell_1,complementarity_iter_ell_inf);
else
    fprintf('H. iters\t CPU Time (m)\t Max. CPU (m)/iter\tMin. CPU (m)/iter\t Comp. res. (ell_1) \tComp. res. (ell_inf)\n');
    fprintf('%d\t\t\t\t%2.2f\t\t%2.2f\t\t\t\t%2.2f\t\t\t\t\t%2.2e\t\t\t\t%2.2e\n',stats.homotopy_iterations,sum(stats.cpu_time)/60,max(stats.cpu_time)/60,min(stats.cpu_time)/60,complementarity_iter_ell_1,complementarity_iter_ell_inf);
end
fprintf('----------------------------------------------------------------------------------------------------------------------\n');
%%
if time_optimal_problem
    T_opt = results.w_opt(model.ind_t_final);
    fprintf('Time optimal problem solved with T_opt: %2.4f.\n',T_opt);
    fprintf('----------------------------------------------------------------------------------------------------------------------\n');
end

%% Output
varargout{1} = results;
varargout{2} = stats;
varargout{3} = model;
varargout{4} = settings;
varargout{5} = solver_initialization;
end