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
function [varargout] = homotopy_solver(varargin)
% homotopy_solver(solver,model,settings,solver_initialization)
solver = varargin{1};
model = varargin{2};
settings = varargin{3};
solver_initialization = varargin{4};

if nargin>4
    model_int = varargin{5};
    settings_int = varargin{6};
end


import casadi.*
%%  unfold data
unfold_struct(settings,'caller')
unfold_struct(solver_initialization,'caller')

comp_res = model.comp_res;
nabla_J_fun = model.nabla_J_fun;
s_elastic_iter = 1;

sigma_k = sigma_0;
x0 = solver_initialization.lbw(1:model.dimensions.n_x);

% lambda00 initialization
if strcmp(settings.pss_mode, 'Stewart')
    lambda00 = full(model.lambda00_fun(x0));
elseif strcmp(settings.pss_mode, 'Step')
%     c_x = full(model.c_fun(x0));
%     lambda00 = [ max(c_x, 0); min(c_x, 0)];
    lambda00 = full(model.lambda00_fun(x0));
end
p_val = [model.p_val(:); lambda00(:)];

% TODO remove try!
try
    complementarity_stats = [full(comp_res(w0, p_val))];
catch
    w0 = w0(1:length(model.w));
    complementarity_stats = [full(comp_res(w0, p_val))];
end
cpu_time = [];
homotopy_iterations = [];
w0_base = w0;
W = [w0];

lbw_h = lbw; ubw_h = ubw;
lbw_h(model.ind_h) = model.h_k(1);
ubw_h(model.ind_h) = model.h_k(1);

%% homotopy loop
complementarity_iter = 1;
ii = 0;

if print_level >= 3
    fprintf('\niter\tsigma\t\tcompl_res\tCPU time\tNLP iters\tstatus\n')
end


while (complementarity_iter) > comp_tol && ii < N_homotopy && sigma_k > sigma_N
    % homotopy parameter update
    if ii == 0
        sigma_k = sigma_0;
    else
        if isequal(homotopy_update_rule,'linear')
            sigma_k = homotopy_update_slope*sigma_k;
        elseif isequal(homotopy_update_rule,'superlinear')
            sigma_k = max(sigma_N,min(homotopy_update_slope*sigma_k,sigma_k^homotopy_update_exponent));
        else
            error('For the homotopy_update_rule please select ''linear'' or ''superlinear''.')
        end  
    end
    p_val(1) = sigma_k;
    if h_fixed_to_free_homotopy
        p_val(3) = 1+(sigma_k*1e4);
    end

    % solve problem with fixed step size
    if h_fixed_iterations && use_fesd  && ii < h_fixed_max_iter
        tic
        results = solver('x0', w0, 'lbx', lbw_h, 'ubx', ubw_h,'lbg', lbg, 'ubg', ubg,'p',p_val);
        cpu_time_iter = toc ;
        w_opt = full(results.x);
        if ~h_fixed_change_sigma
            ii = -1; h_fixed_iterations  = 0;
        end
    else
        tic
        results = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg,'p',p_val);
        cpu_time_iter = toc ;
    end
    if isequal(solver.stats.return_status,'Infeasible_Problem_Detected')
        warning('NLP infeasible: try different mpcc_mode or check problem functions.');
    end

    cpu_time = [cpu_time,cpu_time_iter];
    w_opt = full(results.x);
    w0 = w_opt;
    W = [W,w_opt]; % all homotopy iterations

    % complementarity
    complementarity_iter = full(comp_res(w_opt, p_val));
    complementarity_stats = [complementarity_stats;complementarity_iter];
    % update counter
    ii = ii+1;

    % Verbose
    if print_level >= 3
        fprintf('%d\t\t%2.2e\t%2.2e\t%.3f\t\t%d\t\t\t%s\n',ii, sigma_k, complementarity_iter, ...
            cpu_time_iter, solver.stats.iter_count, solver.stats.return_status);
    end
    %
    %     if complementarity_iter> 1e1 && ii >= ratio_for_homotopy_stop*N_homotopy
    %         error('The homotopy loop is diverging. Try chaning parameters of the homotopy loop or check is the OCP well posed.')
    %         break;
    %     end
end

%% polish homotopy solution with fixed active set.
if polishing_step
    [results] = polishing_homotopy_solution(model,settings,results,sigma_k);
    %     [results] = polishing_homotopy_solution(model,settings,results,sigma_k,solver,solver_initialization);
    complementarity_iter = results.complementarity_iter;
    complementarity_stats = [complementarity_stats;complementarity_iter];
    W = [W,results.w_opt];
end

%% collect stats
results.W = W;
stats.complementarity_stats = complementarity_stats;
stats.cpu_time = cpu_time;
stats.cpu_time_total = sum(cpu_time);
stats.sigma_k = sigma_k;
stats.homotopy_iterations = ii;

%% loop output
varargout{1} = results;
varargout{2} = stats;
varargout{3} = solver_initialization;
end

