clear all
clc
close all

import casadi.*
%% settings
% collocation settings
settings.d = 2;                            % Degree of interpolating polynomial
settings.collocation_scheme = 'radau';     % Collocation scheme: radau or legendre


% control stages


% MPCC settings
settings.solver_name = 'solver_fesd';           % name of solver function.
settings.mpcc_mode = 3;                    % 1 - extact, 2 - smooth  ,3 -relax , 4 - penalty, 5 - elastic mode
settings.s_elastic_max = 1e2;              % upper bound for elastic variables
settings.s_elastic_min = 0;              % upper bound for elastic variables
settings.pointwise_or_integral = 1;        % 1 - pointwise relaxation of the comp., 0 integral relaxation of the cc.
settings.objective_scaling_direct = 1;      % in penalty methods  1: J = J+(1/p)*J_comp (direct)  , 0 : J = p*J+J_comp (inverse)
% Penalty/Relaxation paraemetr
settings.sigma_0 = 1e0;                     % starting smouothing parameter
settings.sigma_N = 1e-10;                   % end smoothing parameter
settings.kappa = 0.1;                      % decrease rate
settings.N_homotopy = 10;
settings.comp_tol = 1e-12;
settings.fesd_complementartiy_mode = 1;
% settings.initial_theta = 1;

%^ IPOPT Settings
opts_ipopt.verbose = false;
opts_ipopt.ipopt.max_iter = 600;
% opts_ipopt.ipopt.max_iter = 100;
opts_ipopt.ipopt.print_level = 5;
opts_ipopt.ipopt.mu_strategy = 'adaptive';
opts_ipopt.ipopt.mu_oracle = 'quality-function';
settings.opts_ipopt = opts_ipopt;
% opts_ipopt.ipopt.tol = 1e-5;
% 1 - scale inverse, 0 scale direct


% finite elements with switch detection
settings.use_fesd = 1;       % turn on moving finite elements algortihm

% step equilibration
settings.gamma_h = 1;                    % how much can h addapt
settings.regularize_h = 1;                 % add quadratic cost term rho_h = (h-hk)^2
settings.rho_h = 0.1;                        % regularization penalty
settings.piecewise_equidistant_grid = 0;
settings.piecewise_equidistant_grid_sigma = 0.1;


%% Time settings

% Tf % Optimal control
settings.time_freezing = 1;
settings.time_optimal_problem = 1;
settings.s_sot_max = 3;

% Time freezing scaling 
settings.use_speed_of_time_variables =  1; % introduce s_tof for e
settings.local_speed_of_time_variable = 1;  

% Grid settings
settings.equidistant_control_grid = 1;
settings.stagewise_clock_constraint = 1;
%% Generate Model
model = car_hystheresis_model();
%% Reformulation of the PSS into a DCS
[model,settings] = model_reformulation_fesd(model,settings);
%% Formulate NLP;
[solver,solver_initalization, model,settings] = create_nlp_fesd(model,settings);

%% Get variables into main workspace
unfold_struct(model,'base');
unfold_struct(settings,'base');
unfold_struct(solver_initalization,'base');
%% Solve NLP
complementarity_stats = [];
cpu_time = [];

% solver_initalization.w0(ind_x(2:n_x:end)) = 49;
[sol,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);
obj = full(sol.f);
w_opt = full(sol.x);
%% Read and plot Result
plot_results_car_hystheresis
