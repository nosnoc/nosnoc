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
settings_bkp = settings;
unfold_struct(settings,'caller');
settings = settings_bkp; % TODO: figure out why unfold settings breaks things.
unfold_struct(model,'caller');
results = extract_results_from_solver(model,settings,results);

% get 00 values
% TODO: Make function for parameters 00
lambda_00 = [];
gamma_00 = [];
p_vt_00 = [];
n_vt_00  = [];
gamma_d00 = [];
delta_d00 = [];
y_gap00 = [];
if settings.dcs_mode == 'CLS'
    % TODO: reconsider this if 0th element has an impulse
    y_gap00 = model.f_c_fun(x0);
    if model.friction_exists
        switch settings.friction_model
          case 'Polyhedral'
            v0 = x0(model.dims.n_q+1:end);
            D_tangent_0 = model.D_tangent_fun(x0);
            v_t0 = D_tangent_0'*v0;
            for ii = 1:model.dims.n_contacts
                ind_temp = model.dims.n_t*ii-(model.dims.n_t-1):model.dims.n_t*ii;
                gamma_d00 = [gamma_d00;norm(v_t0(ind_temp))/model.dims.n_t];
                delta_d00 = [delta_d00;D_tangent_0(:,ind_temp)'*v0+gamma_d00(ii)];
            end
          case 'Conic'
            v0 = x0(model.dims.n_q+1:end);
            v_t0 = model.J_tangent_fun(x0)'*v0;
            for ii = 1:model.dims.n_contacts
                ind_temp = model.dims.n_t*ii-(model.dims.n_t-1):model.dims.n_t*ii;
                v_ti0 = v0(ind_temp);
                gamma_00 = [gamma_00;norm(v_ti0)];
                switch settings.conic_model_switch_handling
                  case 'Plain'
                    % no extra vars
                  case {'Abs','Lp'}
                    p_vt_00 = [p_vt_00;max(v_ti0,0)];
                    n_vt_00 = [n_vt_00;max(-v_ti0,0)];
                end
            end
        end
    end
else
    lambda_00 = full(model.lambda00_fun(x0,model.p_global_val));
end
p_comp = [model.p_val;x0;lambda_00(:);y_gap00(:);gamma_00(:);gamma_d00(:);delta_d00(:);p_vt_00(:);n_vt_00(:)];
complementarity_iter_ell_inf = full(comp_res(results.w_opt,p_comp));
switch dcs_mode
    case 'Step'
    temp = [results.alpha_opt_extended.*results.lambda_n_opt_extended,(1-results.alpha_opt_extended).*results.lambda_p_opt_extended];
    complementarity_iter_ell_1 = sum(temp(:));
    case 'Stewart'
        % TODO: considert cross comps as well in the inf norm
     temp = [results.theta_opt_extended.*results.lam_opt_extended];
    complementarity_iter_ell_1 = sum(temp(:));
    case 'CLS'
%         TODO?
end

%% Verbose
stats.total_time  = total_time;
fprintf('\n');
fprintf('-----------------------------------------------------------------------------------------------\n');
if use_fesd
    fprintf( ['OCP with the FESD ' char(irk_scheme) ' in ' char(irk_representation) ' mode with %d RK-stages, %d finite elements and %d control intervals.\n'],n_s,N_finite_elements(1),N_stages);
else
    fprintf( ['OCP with the Std ' char(irk_scheme) ' in ' char(irk_representation) ' mode with %d RK-stages, %d finite elements and %d control intervals.\n'],n_s,N_finite_elements(1),N_stages);
end

fprintf('---------------------------------------------- Stats summary--------------------------\n');
if sum(stats.cpu_time) < 60
    fprintf('H. iters\t CPU Time (s)\t Max. CPU (s)/iter\tMin. CPU (s)/iter \tComp. res.\n');
    fprintf('%d\t\t\t\t%2.2f\t\t\t%2.2f\t\t\t\t%2.2f\t\t\t\t%2.2e\t\t\t\t%2.2e\n',stats.homotopy_iterations,sum(stats.cpu_time),max(stats.cpu_time),min(stats.cpu_time),complementarity_iter_ell_inf);
else
    fprintf('H. iters\t CPU Time (m)\t Max. CPU (m)/iter\tMin. CPU (m)/iter \tComp. res.\n');
    fprintf('%d\t\t\t\t%2.2f\t\t%2.2f\t\t\t\t%2.2f\t\t\t\t\t%2.2e\t\t\t\t%2.2e \n',stats.homotopy_iterations,sum(stats.cpu_time)/60,max(stats.cpu_time)/60,min(stats.cpu_time)/60,complementarity_iter_ell_inf);
end
fprintf('\n--------------------------------------------------------------------------------------\n');
if time_optimal_problem
    T_opt = results.w_opt(model.ind_t_final);
    fprintf('Time optimal problem solved with T_opt: %2.4f.\n',T_opt);
    fprintf('\n--------------------------------------------------------------------------------------\n');
end

%% Output
varargout{1} = results;
varargout{2} = stats;
varargout{3} = model;
varargout{4} = settings;
varargout{5} = solver_initialization;
end
