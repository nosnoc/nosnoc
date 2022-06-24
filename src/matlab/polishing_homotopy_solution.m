function [results] = polishing_homotopy_solution(model,settings,solver,stats,solver_initalization,results)

import casadi.*
unfold_struct(results,'caller');
unfold_struct(settings,'caller');
unfold_struct(model,'caller');

% x mathaching on points where the algebraic variables are defined
x_modified = x_opt_extended;
x_modified(:,1:n_s+1:end) = [];
% getting the index sets of the algebraic variables
ind_z_extended = reshape(ind_z,n_z,length(ind_z)/n_z);
ind_alpha= [ind_z_extended(1:n_alpha,:)];
ind_lambda0= [ind_z_extended(n_alpha+1:2*n_alpha,:)];
ind_lambda1  = [ind_z_extended(2*n_alpha+1:3*n_alpha,:)];
% evaluating switching futction and cheching the signs
c_eval = [];
for ii = 1:length(x_modified)
    c_eval = [c_eval,full(c_fun(x_modified(:,ii)))];
end
eps_sigma = stats.sigma_k*10;
ind_negative = c_eval < eps_sigma;
ind_positive = c_eval > eps_sigma ;
ind_sliding = abs(c_eval) <= eps_sigma ;
% c > 0 , lambda0 = 0, alpha = 1;
ind_lambda0_fixed = ind_lambda0(ind_positive);
ind_alpha1_fixed = ind_alpha(ind_positive);
% c < 0 , lambda1 = 0, alpha = 0;
ind_lambda1_fixed = ind_lambda1(ind_negative);
ind_alpha0_fixed = ind_alpha(ind_negative);
% c = 0 , lambda1 = 0, lambda0 = 0;
% her we must make clear is it sliding or is it corssing
ind_lambda1_sliding = ind_lambda1(ind_sliding);
ind_lambda0_sliding = ind_lambda0(ind_sliding);
% set appropiate bounds to zero
unfold_struct(solver_initalization,'caller')
lbw(ind_lambda1_fixed(:)) = 0;  ubw(ind_lambda1_fixed(:)) = 0;
lbw(ind_alpha0_fixed(:)) = 0;  ubw(ind_alpha0_fixed(:)) = 0;

lbw(ind_lambda0_fixed(:)) = 0;  ubw(ind_lambda0_fixed(:)) = 0;
lbw(ind_alpha1_fixed(:)) = 1;  ubw(ind_alpha1_fixed(:)) = 1;
%
lbw(ind_lambda1_sliding(:)) = 0;  ubw(ind_lambda1_sliding(:)) = 0;
lbw(ind_lambda0_sliding(:)) = 0;  ubw(ind_lambda0_sliding(:)) = 0;
% solve nlp
results = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg,'p',0);
% extract all results from w_opt (create a function for this, that maps them to the results
results = extract_results_from_solver(model,settings,results);
end

