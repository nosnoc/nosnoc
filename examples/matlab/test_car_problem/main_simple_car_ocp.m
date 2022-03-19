clear all
clc
close all

import casadi.*

%% settings
[settings] = default_settings_fesd();  %% Optionally call this function to have an overview of all options. 
                                        % Missing settings are anyway filled in latter with their respecitve values.

% collocation settings
settings.d = 2;                            % Degree of interpolating polynomial
settings.collocation_scheme = 'radau';     % Collocation scheme: radau or legendre
settings.collocation_scheme = 'legendre';     % Collocation scheme: radau or legendre

% MPCC settings
settings.solver_name = 'solver_fesd';           % name of solver function.
settings.mpcc_mode = 8;                    % 1 - extact, 2 - smooth  ,3 -relax , 4 - penalty, 5 - elastic mode
settings.s_elastic_max = 1e2;              % upper bound for elastic variables
settings.s_elastic_min = 0;              % upper bound for elastic variables
settings.s_elastic_0 = 1;              % upper bound for elastic variables
settings.pointwise_or_integral = 1;        % 1 - pointwise relaxation of the comp., 0 integral relaxation of the cc.
settings.objective_scaling_direct = 1;      % in penalty methods  1: J = J+(1/p)*J_comp (direct)  , 0 : J = p*J+J_comp (inverse)
settings.lp_initalization = 0;
settings.initial_theta = 0.5*0;
settings.initial_lambda = 0.5*0;
settings.initial_mu = 0.5*0;
% Penalty/Relaxation paraemetr
settings.sigma_0 = 1e2;                     % starting smouothing parameter
settings.sigma_N = 1e-14;                   % end smoothing parameter
settings.kappa = 0.1;                      % decrease rate
settings.N_homotopy = ceil(abs(log(settings.sigma_N/settings.sigma_0)/log(settings.kappa)))+1 ;% number of steps

settings.N_homotopy = 8;

% settings.N_homotopy = 1;% number of steps
settings.comp_tol = 1e-14;
settings.cross_complementarity_mode = 3;
% settings.initial_theta = 1;

%^ IPOPT Settings
opts_ipopt.verbose = false;
opts_ipopt.ipopt.max_iter = 600;
opts_ipopt.ipopt.max_iter = 1000;
opts_ipopt.ipopt.print_level = 5;
opts_ipopt.ipopt.mu_strategy = 'adaptive';
opts_ipopt.ipopt.mu_oracle = 'quality-function';
settings.opts_ipopt = opts_ipopt;


% finite elements with switch detection
settings.use_fesd = 1;       % turn on moving finite elements algortihm
settings.store_all_homotopy_iterates = 1;
% step equilibration
settings.gamma_h = 1;                    % how much can h addapt
% regularize_h --> heuristic_step_equilibration 
settings.regularize_h = 1;                 % add quadratic cost term rho_h = (h-hk)^2
settings.rho_h = 1;                        % regularization penalty
% piecewise_equidistant_grid -- > step_equilibration
% step_equilibration_mode = 1; 1...6

settings.heuristic_step_equilibration = 0;
settings.step_equilibration = 1;
settings.step_equilibration_mode = 1;
settings.step_equilibration_sigma = 0.01;
settings.treat_step_equilibration_with_mpcc_method = 0;
settings.step_equilibration_penalty = 1e2;

settings.convex_sigma_rho_constraint = 1;
settings.sigma_penalty = 1e2;
settings.rho_penalty = 1;

%% Time settings

% Tf % Optimal control
settings.time_freezing = 0;
settings.time_optimal_problem = 1;

% Time freezing scaling 
settings.use_speed_of_time_variables =  1; % introduce s_tof for e
settings.local_speed_of_time_variable = 1;  

% Grid settings
settings.equidistant_control_grid = 1;
settings.stagewise_clock_constraint = 0;

[results,stats,solver_initalization,settings,model] = solve_simple_car_ocp(settings);

%% Read and plot Result
if results.status == 1
fprintf(results.messages.solver_message)
fprintf(results.messages.time_message)
fprintf(results.messages.objective_message)
plot_results_simple_car(results,settings,model,stats)
else
    fprintf('diverged.\n')
end


model.n_cross_comp
%%
figure
stairs(results.w_opt(model.ind_h))
grid on

