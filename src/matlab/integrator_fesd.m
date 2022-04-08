%
%    This file is part of NOS-NOC.
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
function [varargout] = integrator_fesd(model,settings)
import casadi.*

%%  unfold data
% use_previous_solution_as_initial_guess = 1;

unfold_struct(settings,'caller')
unfold_struct(model,'caller')

results.t_grid = [];
results.x_res = [];
stats = [];

%% Settings of integration
[model] = refine_model_integrator(model,settings);
% [settings] = refine_settings_integrator(settings,model)

%% Create solver functions for integrator step

[solver,solver_initalization, model,settings] = create_nlp_fesd(model,settings);
unfold_struct(settings,'caller')
unfold_struct(model,'caller')
unfold_struct(solver_initalization,'caller')
%% Initalization
x_res = [x0];
theta_res = [];
lambda_res = [];
mu_res = [];

% Output including all RK-stage values
x_res_extended = [x0];
lambda_res_extended = [];
theta_res_extended = [];
mu_res_extended = [];

% if c_{n_s} \neq 1 output also the boundary lambda and mu
theta_boundary_res = [];
lambda_boundary_res = [];
mu_boundary_res = [];
% step size
h_vec = [];
% statse
complementarity_stats  = [];
homotopy_iteration_stats = [];
time_per_iter = [];
%% Main simulation loop
for ii = 1:N_sim
    [sol,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);
    time_per_iter = [time_per_iter; stats.cpu_time_total];
    % verbose
    if stats.complementarity_stats(end) > 1e-3
        error('NLP Solver did not converge for the current FESD problem. \n')
    end
    if print_level >=2
        fprintf('Integration step %d / %2.0f (%2.3f s / %2.3f s) converged in %2.3f seconds. \n',ii,N_sim,ii*T,N_sim*T,time_per_iter(end));
    end
    % Store differentail states
    w_opt = full(sol.x);
    diff_states = w_opt(ind_x);
    alg_states = w_opt(ind_z);
    alg_states_boundary = w_opt(ind_boundary);

    if ~right_boundary_point_explicit && use_fesd
        try
            alg_states_boundary = reshape(alg_states_boundary,(n_z-n_theta),N_finite_elements(1)-1);
        catch

        end
    end
    % step-size
    h_opt = w_opt(ind_h); 
    % differential
    x_opt_extended = w_opt(ind_x);
    x_opt_extended  = reshape(x_opt_extended,n_x,length(x_opt_extended)/n_x);
    % only bounadry value
    if isequal(irk_representation,'integral') || lift_irk_differential
        x_opt  = x_opt_extended(:,1:n_s+1:end);
    end
    %algebraic
    alg_states_extended = reshape(alg_states,n_z,length(alg_states)/n_z);
    theta_opt_extended = alg_states_extended(1:n_theta,:);
    lambda_opt_extended=  alg_states_extended(n_theta+1:2*n_theta,:);
    mu_opt_extended  = alg_states_extended(end-n_simplex+1:end,:);
    
    % Only boundary value
    theta_opt  = theta_opt_extended(:,1:n_s:end);
    lambda_opt  = lambda_opt_extended(:,1:n_s:end);
    mu_opt  = mu_opt_extended(:,1:n_s:end);

    % update initial guess and inital value
    x0 = x_opt(:,end);
    solver_initalization.lbw(1:n_x) = x0;
    solver_initalization.ubw(1:n_x) = x0;
    if use_previous_solution_as_initial_guess
        solver_initalization.w0 = w_opt;
    end

    % Store data
    % step-size
    if use_fesd
        h_vec = [h_vec;h_opt];
    else
        h_vec = [h_vec;h*ones(N_stages,1)];
    end
    %differntial
    x_res = [x_res, x_opt(:,end-N_finite_elements(1)*N_stages+1:end)];
    x_res_extended = [x_res_extended,x_opt_extended(:,2:end)];
    % algebraic
    theta_res = [theta_res, theta_opt];
    lambda_res = [lambda_res, lambda_opt];
    mu_res = [mu_res, mu_opt];
    theta_res_extended = [theta_res_extended,theta_opt_extended ];
    lambda_res_extended = [lambda_res_extended,lambda_opt_extended];
    mu_res_extended = [mu_res_extended,mu_opt_extended];
    if ~right_boundary_point_explicit && use_fesd
        lambda_boundary_res = [lambda_boundary_res,alg_states_boundary(1:n_theta,:)];
        mu_boundary_res = [mu_boundary_res,alg_states_boundary(end-n_simplex+1:end,:) ];
    end
    %stats
    complementarity_stats  = [complementarity_stats; stats.complementarity_stats(end)];
    homotopy_iteration_stats = [homotopy_iteration_stats;stats.homotopy_iterations];

end
total_time = sum(time_per_iter);
%% Verbose
fprintf('\n');
fprintf('-----------------------------------------------------------------------------------------------\n');
fprintf( ['Integration procedure with the FESD ' irk_scheme ' with %d stages completed in %2.3f seconds.\n'],n_s,total_time);
fprintf( ['IRK representation: ' irk_representation '.\n']);
fprintf('Total integration steps: %d, with nominal step-size h =  %2.3f and %d finite elements.\n',N_sim,h_sim,N_stages);
if additional_residual_ingeration_step
    fprintf('--> + additional residual step to reach T_sim with  T_residual =  %2.3f.\n',T_residual);
end
fprintf('-------------------------------------------------------------------------------------\n');
fprintf('Total CPU time: %2.3f s, N_stg =%d, N_FE = %d.\n',sum(time_per_iter),N_stages,N_finite_elements(1));
fprintf('Max iteration time: %2.3f s, min iteration time: %2.3f s.\n',max(time_per_iter),min(time_per_iter));
fprintf('Max complementarity residual: %2.3e, min complementarity residual: %2.3e.\n',max(complementarity_stats),min(complementarity_stats));
fprintf('-----------------------------------------------------------------------------------------------\n\n');
%% Output

results.h_vec  = h_vec;
results.t_grid = cumsum([0;h_vec])';

results.x_res  = x_res;
results.theta_res = theta_res;
results.lambda_res = lambda_res;
results.mu_res = mu_res;

% Output all stage values as well.
results.x_res_extended  = x_res_extended;
results.theta_res_extended  = theta_res_extended;
results.lambda_res_extended  = lambda_res_extended;
results.mu_res_extended  = mu_res_extended;
% Output boundary points
results.lambda_boundary_res= lambda_boundary_res;
results.mu_boundary_res= mu_boundary_res;

stats.complementarity_stats   = complementarity_stats;
stats.time_per_iter = time_per_iter;
stats.homotopy_iteration_stats = homotopy_iteration_stats;

varargout{1} = results;
varargout{2} = stats;
varargout{3} = model;
varargout{4} = settings;
end

