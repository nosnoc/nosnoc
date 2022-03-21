function [results,stats] = integrator_fesd(model,settings)
import casadi.*

%%  unfold data
% use_previous_solution_as_initial_guess = 1;

unfold_struct(settings,'caller')
unfold_struct(model,'caller')

results.t_grid = [];
results.x_res = [];

stats = [];

%% Settings

if N_finite_elements > 1
    warning(' N_finite_elements > 1, setting to 1 in simulation problem.')
    settings.N_finite_elements = 1;
end

if N_stages > 3
    warning(' N_stages > 3, smaller values might lead to better perfomance.')
end

if N_stages*h ~= T
    model.T = N_stages*h;
end



%% Create solver functions for integrator step
settings.opts_ipopt.ipopt.print_level=0;
[solver,solver_initalization, model,settings] = create_nlp_fesd(model,settings);
unfold_struct(settings,'caller')
unfold_struct(model,'caller')
unfold_struct(solver_initalization,'caller')
% TODO: figure out how to deal with the last T if T_sim/N_sim is not a whole number;
%% Initalization
x_res = [x0];
theta_res = [];
lambda_res = [];
mu_res = [];
h_vec = [];
complementarity_stats  = [];
homotopy_iteration_stats = [];
if exist('N_sim')
    %     || isempty(N_sim)
    if N_sim*(N_stage*h) ~= T_sim
        warning(' The given N_sim does not result in the desired, given T_sim. Do we fix it automatically now?')
    end
else
    N_sim = round(T_sim/(h*N_stages));
end
h_vec = [];
time_per_iter = [];
%% Main simulation loop
for ii = 1:N_sim
    tic
    [sol,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);
    time_per_iter = [time_per_iter; toc];
%     if stats.complementarity_stats(end) > 1e-3
%         error('NLP Solver did not converge for the current FESD Problem. \n')
%     end
    fprintf('Integration step %d of %d (%2.2.f s /%2.2.f s )converged in %2.2f seconds. \n',ii,N_sim,ii*T,T_sim,time_per_iter(end));

    % Store differentail states
    w_opt = full(sol.x);
    diff_states = w_opt(ind_x);
    alg_states = w_opt(ind_z);
    h_opt = w_opt(ind_h);
    % differential and algebraic states
    x_opt = [];
    theta_opt = [];
    lambda_opt = [];
    mu_opt = [];

    for i = 1:n_x
        eval( ['x' num2str(i) '_opt = diff_states(' num2str(i) ':n_x+n_x*d:end);']);
        eval(['x_opt = [x_opt,; transpose(x' num2str(i) '_opt)];'])
    end

 % convex multiplers
    for i = 1:n_theta
        eval( ['theta' num2str(i) '_opt = alg_states(' num2str(i) ':n_z+n_z*(d-1):end);']);
        eval(['theta_opt = [theta_opt,; transpose(theta' num2str(i) '_opt)];'])
    end
% lambdas
    for i = 1:n_theta
        eval( ['lambda' num2str(i) '_opt = alg_states(' num2str(i+n_theta) ':n_z+n_z*(d-1):end);']);
        eval(['lambda_opt = [lambda_opt,; transpose(lambda' num2str(i) '_opt)];'])
    end
% mu
    for i = 1:n_simplex   
        eval(['mu' num2str(i) '_opt = alg_states(' num2str(i+2*n_theta) ':n_z+n_z*(d-1):end);']);
        eval(['mu_opt = [mu_opt,; transpose(mu' num2str(i) '_opt)];'])
    end


    x0 = x_opt(:,end);
    solver_initalization.lbw(1:n_x) = x0;
    solver_initalization.ubw(1:n_x) = x0;

    if use_previous_solution_as_initial_guess
        solver_initalization.w0 = w_opt;
    end
    %% Store data
    h_vec = [h_vec;h_opt];
    complementarity_stats  = [complementarity_stats; stats.complementarity_stats(end)];
    homotopy_iteration_stats = [homotopy_iteration_stats;stats.homotopy_iterations];
    x_res = [x_res, x_opt(:,end-N_stages+1:end)];
    theta_res = [theta_res, theta_opt(:,end-N_stages+1:end)];
    lambda_res = [lambda_res, lambda_opt(:,end-N_stages+1:end)];
    mu_res = [mu_res, mu_opt(:,end-N_stages+1:end)];
end
total_time = sum(time_per_iter);
%% Verbose
fprintf('\n');
fprintf('-------------------------------------------------------------------------------------\n');
switch collocation_scheme
    case 'radau'
        fprintf('Integration procedure with the FESD Radau IIA - %d scheme completed in %2.3f seconds.\n',2*d-1,total_time);
    case 'legendre'
        fprintf('Integration procedure with the FESD Gauss-Legendre - %d scheme completed in %2.3f seconds.\n',2*d,total_time);
end
fprintf('Total integration steps: %d, with nominal step-size h =  %2.3f and %d stages.\n',N_sim,h,N_stages);
fprintf('-------------------------------------------------------------------------------------\n');
fprintf('Max iteration time: %2.3f seconds, min iteration time: %2.3f seconds.\n',max(time_per_iter),min(time_per_iter));
fprintf('Max complementarity residual: %2.3e, min complementarity residual: %2.3e.\n',max(complementarity_stats),min(complementarity_stats));
fprintf('-------------------------------------------------------------------------------------\n');
%% Output

results.t_grid = cumsum([0;h_vec])';
results.x_res  = x_res;
results.theta_res = theta_res;
results.lambda_res = lambda_res;
results.mu_res = mu_res;
results.h_vec  = h_vec;
stats.complementarity_stats   = complementarity_stats;
stats.time_per_iter = time_per_iter;
stats.homotopy_iteration_stats = homotopy_iteration_stats;
end

