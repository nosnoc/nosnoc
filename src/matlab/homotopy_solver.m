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
% homotopy_solver(solver,model,settings,solver_initalization)
solver = varargin{1};
model = varargin{2};
settings = varargin{3};
solver_initalization = varargin{4};

if nargin>4
    model_int = varargin{5};
    settings_int = varargin{6};
end


import casadi.*
%%  unfold data
unfold_struct(settings,'caller')
unfold_struct(solver_initalization,'caller')
p_val = model.p_val;

comp_res = model.comp_res;
nabla_J_fun = model.nabla_J_fun;
s_elastic_iter = 1;

sigma_k = sigma_0;
x0 = model.x0;

try
    complementarity_stats = [full(comp_res(w0))];
catch
    w0 = w0(1:length(model.w));
    complementarity_stats = [full(comp_res(w0))];
end
cpu_time = [];
homotopy_iterations = [];
w0_base = w0;
W = [w0];

lbw_h = lbw; ubw_h = ubw;
lbw_h(model.ind_h) = model.h_k(1);
ubw_h(model.ind_h) = model.h_k(1);

%% homtopy loop
complementarity_iter = 1;
ii = 0;
vf_resiudal = 0;




while (complementarity_iter+vf_resiudal) > comp_tol && ii < N_homotopy
    % homotopy parameter update
    if ii == 0
        sigma_k = sigma_0;
    else
        sigma_k = kappa*sigma_k;
    end
    p_val(1) = sigma_k;
    if h_fixed_to_free_homotopy
        p_val(3) = 1+(sigma_k*1e4);
    end
    %     if integrator_forward_sweep_procedure
    %         % do not do it in first iteration if kinematic presolve was done
    %         if (virtual_forces_kinematic_iteration && ii >0) && ii >=0
    %           u_sim = w0(model.ind_u);
    %           model_int.T_sim = model.h;
    %           settings_int.mpcc_mode = 4;
    %           settings_int.print_level = 2;
    %           u_sim = reshape(u_sim,model.n_u,model.N_stages);
    %           [results,stats] = integrator_fesd(model_int,settings_int,u_sim);
    % %           w0 = integrator_forward_sweep(model_int,solver_int,solve_initalization_int);
    %         end
    %     end
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

    cpu_time = [cpu_time,cpu_time_iter];
    w_opt = full(results.x);
    w0 = w_opt;
    W = [W,w_opt]; % all homotopy iterations

    % complementarity
    complementarity_iter = full(comp_res(w_opt));
    complementarity_stats = [complementarity_stats;complementarity_iter];
    if virtual_forces
        vf_resiudal = full(model.J_virtual_froces_fun(w_opt));
    end
    % update counter
    ii = ii+1;

    % Verbose
    if print_level>=3
        fprintf('-----------------------------------------------------------------------------------------------\n');
        fprintf('Homotopy iteration : %d / %d, with sigma = %2.2e completed.\n',ii,N_homotopy,sigma_k);
        fprintf('Complementarity resiudal: %2.2e.\n',complementarity_iter);
        if model.n_u >0
            fprintf('CPU time of iteration: %2.2f s.\n',cpu_time_iter);
            fprintf('Objective function value: %2.4e.\n',cpu_time_iter);
            if time_optimal_problem
                fprintf('Final time T_opt: %2.4f.\n',w_opt(model.ind_t_final));
            end
            if virtual_forces
                fprintf('Virtual forces residual: %2.2e.\n',vf_resiudal);
            end
        else
            fprintf('CPU time of iteration: %2.2f s.\n',cpu_time_iter);
        end
        fprintf('-----------------------------------------------------------------------------------------------\n');
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
    %     [results] = polishing_homotopy_solution(model,settings,results,sigma_k,solver,solver_initalization);
    complementarity_iter = results.complementarity_iter;
    complementarity_stats = [complementarity_stats;complementarity_iter];
    W = [W,results.w_opt];
end

%% collcet stats
results.W = W;
stats.complementarity_stats = complementarity_stats;
stats.cpu_time = cpu_time;
stats.cpu_time_total = sum(cpu_time);
stats.sigma_k = sigma_k;
stats.homotopy_iterations = ii;

%% loop output
varargout{1} = results;
varargout{2} = stats;
varargout{3} = solver_initalization;
end

