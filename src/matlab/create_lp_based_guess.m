function [theta_guess,lambda_guess,mu_guess] = create_lp_based_guess(model)
% creates a guess for the algebraic variables by solving the Stewart LP for
% x = x0
import casadi.*

unfold_struct(model,'caller');

%% Guss
initial_theta = 1;
theta_guess = initial_theta*ones(n_theta,1)*1;
lambda_guess = theta_guess;
mu_guess = initial_theta*ones(n_simplex,1);
%% objective gradient of the LP
try
    g_guess = full(g_Stewart_fun(x0));
catch
    g_guess = full(g_Stewart_fun(x0,u0));
end
%% equality constraint for LP
Aeq = (E');
beq = ones(n_simplex,1);


 %% Ipopt setting for efficient LP solving
    opts_ipopt_lp.verbose = false;
    opts_ipopt_lp.ipopt.max_iter = 200;
    opts_ipopt_lp.ipopt.fixed_variable_treatment = 'make_constraint';
    opts_ipopt_lp.ipopt.mehrotra_algorithm = 'yes';                            % barrier update strategy that works very well for convex QPs
    opts_ipopt_lp.ipopt.jac_c_constant = 'yes';
    opts_ipopt_lp.ipopt.jac_d_constant = 'yes';
    opts_ipopt_lp.ipopt.hessian_constant = 'yes';
    opts_ipopt_lp.ipopt.honor_original_bounds = 'yes';
    tol_ipopt = 1e-16;
    opts_ipopt_lp.ipopt.tol = tol_ipopt;
    opts_ipopt_lp.ipopt.dual_inf_tol = tol_ipopt;
    opts_ipopt_lp.ipopt.dual_inf_tol = tol_ipopt;
    opts_ipopt_lp.ipopt.compl_inf_tol = tol_ipopt;
    opts_ipopt_lp.ipopt.bound_relax_factor = 0.0;
    opts_ipopt_lp.ipopt.print_level = 5;
    
    %% Create LP 
    theta_lp = MX.sym('theta_lp',n_theta);

    f_lp = g_guess'*theta_lp;
    g_lp = Aeq*theta_lp - beq;
    lbw_lp = 0*ones(n_theta,1);
    ubw_lp = inf*ones(n_theta,1);
    prob_lp = struct('f', f_lp, 'x', theta_lp, 'g',g_lp);
    solver_lp = nlpsol('slover_lp', 'ipopt', prob_lp,opts_ipopt_lp);
    % Solve LP 
    sol_lp = solver_lp('x0', theta_guess, 'lbx', lbw_lp, 'ubx',ubw_lp,'lbg', beq*0, 'ubg', beq*0);
    % Read solution
    theta_guess = full(sol_lp.x);
    lambda_guess = -full(sol_lp.lam_x);
    mu_guess = full(sol_lp.lam_g);
    %% Quadprog sol
    %     %     LB = zeros(n_theta,1);
%     %     [theta_opt,~,~,~,multiplers]  = linprog(h_guess ,[],[],Aeq,beq,LB,[],[]);
%     %     %     theta_guess = theta_guess;
%     %         lambda_guess = multiplers.lower;
%     %     %     mu_guess = multiplers.eqlin;

end

