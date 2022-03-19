
clear all
clc
close all

import casadi.*

% digits(64)

%% Benchmark settings
d_vec = [2];
% M_vec = [50;100;200;500;1000;1500;2000;2500;3000;4000];

% M_vec = [50;100;150;200;];
% M_vec = [50;100;200;500;1000;1500;2000;2500;3000;4000;5000;10000];
M_vec  = 200;
% d_vec = 1;
% M_vec = 500;
% M_vec = [100;200;500;1000;1500;2000];


legend_str ={'Implicit Euler','IRK Radau 3','IRK Radau 5','IRK Radau 7','IRK Radau 9'};
% legend_str ={'IRK Radau 3','IRK Radau 5'};
legend_str = [legend_str(d_vec)];

%% settings
% collocation settings
settings.d = 2;                            % Degree of interpolating polynomial
settings.collocation_scheme = 'radau';     % Collocation scheme: radau or legendre
settings.initial_guess = 1;              %  1  - zeros, 2 - time interpolation , 3- sinus/cos

% MPCC settings
settings.solver_name = 'solver_mfe';           % name of solver function.
settings.mpcc_mode = 4;                    % 1 - extact, 2 - smooth  ,3 -relax , 4 - penalty, 5 - l1 penalty
settings.pointwise_or_integral = 0;        % 1 - pointwise relaxation of the comp., 0 integral relaxation of the cc.
settings.integral_full = 1;                % 1 - all_complementarites_of the NLP into single constraint
settings.l1_scaling = 0;
% Penalty/Relaxation paraemetr
settings.sigma_0 = 1e-1;                     % starting smouothing parameter
settings.sigma_N = 1e-4;                   % end smoothing parameter
settings.kappa = 0.1;                      % decrease rate
settings.N_homotopy = ceil(abs(log(settings.sigma_N/settings.sigma_0)/log(settings.kappa)))+1 ;% number of steps


settings.comp_tol = 1e-12;
settings.sigma_0 = 1e0;                     % starting smouothing parameter
settings.N_homotopy = 30 ;% number of steps
settings.kappa = 0.5;                      % decrease rate

%^ IPOPT Settings
opts_ipopt.verbose = false;
opts_ipopt.ipopt.max_iter = 1500;
opts_ipopt.ipopt.mu_strategy = 'adaptive';
opts_ipopt.ipopt.mu_oracle = 'quality-function';
opts_ipopt.ipopt.tol = 1e-16;
opts_ipopt.ipopt.dual_inf_tol = 1e-16;              % Desired threshold for the dual infeasibility.
opts_ipopt.ipopt.constr_viol_tol = 1e-16;           % Desired threshold for the constraint violation.
opts_ipopt.ipopt.compl_inf_tol = 1e-16;           % Desired threshold for the complementarity conditions.

settings.opts_ipopt = opts_ipopt;



% moving finite elements
settings.moving_finite_elements = 1;       % turn on moving finite elements algortihm
settings.gamma_h = 1-1e-16;                    % how much can h addapt
% settings.gamma_h = 0.8;                    % how much can h addapt
settings.regularize_h = 0;                 % add quadratic cost term rho_h = (h-hk)^2
settings.rho_h = 1;                        % regularization penalty


%% Time settings

omega = 2*pi;
x_star = [exp(1);0];
s_star = [exp(2)  0; exp(2)*2*omega exp(2)];


T = 2;                            % time budget of transformed pseudo time
T_sim = T;



model.N_stages = N;

%% Generate Model

% model = bouncing_ball(model);

%% Reformulation of the PSS into a DCS
% [model,settings] = model_reformulation_mfe(model,settings);

%% Formulate NLP;
% [solver,solver_initalization, model,settings] = create_nlp_mfe(model,settings);

% %% Get variables into main workspace
% unfold_struct(model,'base');
% unfold_struct(settings,'base');
% unfold_struct(solver_initalization,'base');

%% Solve NLP
errors_all_experiments = [];
errors_switch_detection_1_all_experiments = [];
errors_switch_detection_2_all_experiments = [];
complementarity_all_experiments = [];
nominal_h_all_experiments = [];

% optimal siwtch points
t1_star =1;




for i = 1:length(d_vec)
    d = d_vec(i);
    n_col = 2*(d+1); % number of collocation points per 2 finite elements
    
    settings.d = d; % update collocation order
    
    % store data for fixed d and variable M/h
    errors_current_experiment = [];
    errors_switch_detection_1_current_experiment = [];
    complementarity_current_experiment = [];
    nominal_h_current_experiment = [];
    
    for  j = 1:length(M_vec)
        M = M_vec(j);
        N_sim = round(M/n_col); % total number of simulation intevals
        h = T_sim/(2*N_sim);    % nominal step lenght of a single finite elemnt;
        
        % update step size
        model.T = 2*h;
        model.h = h;
        
        % generate new model with updated settings;
        model = oscilator(model);
        [model,settings] = model_reformulation_mfe(model,settings);
        [solver,solver_initalization, model,settings] = create_nlp_mfe(model,settings); % generate new NLP with updated settings
        
        % Update Data
        unfold_struct(model,'base');
        unfold_struct(settings,'base');
        unfold_struct(solver_initalization,'base');
        
        
        
        sigma_k = sigma_0;
        x0 = model.x0;
        
        complementarity_stats = [];
        cpu_time = [];
        diff_states = [x0];
        controls = [];
        homotopy_iteration_stats = [];
        w0_base = w0;
        
        for jj  = 1:N_sim
            % simulation for current experimetn
            tic
            solver_initalization.lbw(1:n_x) = x0;
            solver_initalization.ubw(1:n_x) = x0;
            %             w0 = w0_base;
            [sol,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);
            homotopy_iteration_stats = [homotopy_iteration_stats;stats.homotopy_iterations];
            w_opt = full(sol.x);
            solver_initalization.w0 = w_opt;
            % read result
            x_opt = w_opt(ind_x);
            x_opt = reshape(x_opt,n_x,N_stages*(d+1)+1);
            
            u_opt = w_opt(ind_u);
            
            diff_states = [diff_states,x_opt(:,[(d+1)+1,N_stages*(d+1)+1])];
            controls = [controls,u_opt(1:2),u_opt(3:end)];
            x0 = diff_states(:,end);
            
            f_obj = full(sol.f);
            
            
            complementarity_iter = full(comp_res(w_opt));
            complementarity_stats = [complementarity_stats;complementarity_iter ];
        end
        % stats for current step size
        h_opt = controls(2,:)';
        tgrid = (cumsum([0;h_opt]));
        tgrid_z = cumsum(h_opt);
        
      
        
        % numerical error
        x_mfe = diff_states(1:2,end);
        error_x = norm(x_mfe-x_star,inf);
        t1_mfe = sum(h_opt(1:N_sim));
  
        error_t1 = abs(1-t1_mfe);
        % complementarity
        max_complementarity_exp = max(complementarity_stats);
        
        % save date current experiemnt
        errors_current_experiment = [errors_current_experiment,error_x];
        fprintf('Error with (h = %2.5f, M = %d, d = %d ) is %5.2e : \n',h,M,d,error_x);
        fprintf('Complementarity residual %5.2e : \n',max_complementarity_exp);
        
        complementarity_current_experiment = [complementarity_current_experiment,max_complementarity_exp];
        nominal_h_current_experiment = [nominal_h_current_experiment,h];
        
        errors_switch_detection_1_current_experiment = [errors_switch_detection_1_current_experiment,error_t1];
    end
    % all errors (the vectors of every scheme) an n(d_vec) x n(M_vec) matrices
    
    % numerical error
    errors_all_experiments = [errors_all_experiments;errors_current_experiment];
    
    % switch detection
    errors_switch_detection_1_all_experiments = [errors_switch_detection_1_all_experiments;errors_switch_detection_1_current_experiment];
    
    nominal_h_all_experiments = [nominal_h_all_experiments;nominal_h_current_experiment];
    complementarity_all_experiments = [complementarity_all_experiments;complementarity_current_experiment];
end

%%
diff_states = w_opt(ind_x);
controls = w_opt(ind_u);
alg_states = w_opt(ind_z);

% differential states
for i = 1:n_x
    eval( ['x' num2str(i) '_opt = diff_states(' num2str(i) ':n_x+n_x*d:end);']);
end

%
figure
plot(x1_opt,x2_opt);