clear all
clc
close all

import casadi.*
%% Settings
% collocation settings
settings.d = 4;                            % Degree of interpolating polynomial
settings.d = 3;                            % Degree of interpolating polynomial
settings.collocation_scheme = 'radau';     % Collocation scheme: radau or legendre
settings.lp_initalization = 0;

% control stages

% MPCC settings
settings.solver_name = 'solver_fesd';           % name of solver function.
settings.mpcc_mode = 5;                    % 1 - extact, 2 - smooth  ,3 -relax , 4 - penalty, 5 - elastic mode
settings.s_elastic_max = 1e1;              % upper bound for elastic variables
settings.pointwise_or_integral = 1;        % 1 - pointwise relaxation of the comp., 0 integral relaxation of the cc.
settings.objective_scaling_direct = 1;                   % in penalty methods  1: J = J+(1/p)*J_comp (direct)  , 0 : J = p*J+J_comp (inverse)
% Penalty/Relaxation paraemetr
settings.sigma_0 = 1e0;                     % starting smouothing parameter
settings.sigma_N = 1e-5;                   % end smoothing parameter
settings.sigma_N = 1e-12;                   % end smoothing parameter
settings.kappa = 0.1;                      % decrease rate
settings.N_homotopy = ceil(abs(log(settings.sigma_N/settings.sigma_0)/log(settings.kappa)))+1 ;% number of steps
settings.comp_tol = 1e-15;

%^ IPOPT Settings
opts_ipopt.verbose = false;
opts_ipopt.ipopt.max_iter = 1500;
% opts_ipopt.ipopt.max_iter = 100;
opts_ipopt.ipopt.print_level = 0;
opts_ipopt.ipopt.print_level = 5;
opts_ipopt.ipopt.mu_strategy = 'adaptive';
opts_ipopt.ipopt.mu_oracle = 'quality-function';
settings.opts_ipopt = opts_ipopt;

% FESD
settings.use_fesd = 1;       
settings.moving_finite_elements = 1;       
settings.gamma_h = 1;                    % how much can h addapt
settings.regularize_h = 1;                 % add quadratic cost term rho_h = (h-hk)^2
settings.rho_h = 0.1;                        % regularization penalty
settings.piecewise_equidistant_grid = 0;
settings.fesd_complementartiy_mode = 1;
settings.delta_h_regularization = 0;
settings.piecewise_equidistant_grid_sigma = 0.1;
settings.equidistant_control_grid = 1;

%% Time settings

% time freezing scaling 
settings.use_speed_of_time_variables =  1; % introduce s_tof for e
settings.local_speed_of_time_variable = 1;  
% Grid settings
settings.equidistant_control_grid = 1;

settings.time_freezing = 1;
settings.stagewise_clock_constraint = 1;
% settings.impose_terminal_phyisical_time = 1;

%% Generate Model
model = throwing_ball_model();
%% Reformulation of the PSS into a DCS
[model,settings] = model_reformulation_fesd(model,settings);
%% Formulate NLP;
if 1
    [solver,solver_initalization, model,settings] = create_nlp_fesd(model,settings);
else
    settings.store_all_homotopy_iterates = 1;
    settings.moving_finite_elements = 1;
    settings.l1_scaling = 0;
    [solver,solver_initalization, model,settings] = create_nlp_mfe_develop(model,settings);
end
%% Get variables into main workspace
unfold_struct(model,'base');
unfold_struct(settings,'base');
unfold_struct(solver_initalization,'base');
%% Solve NLP
[sol,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);
obj = full(sol.f);
w_opt = full(sol.x);
%% Read and plot Result
plot_results_thowring_ball
[x1_opt(1),x2_opt(1),x3_opt(1),x4_opt(1),x5_opt(1)]
[x1_opt(end),x2_opt(end),x3_opt(end),x4_opt(end),x5_opt(end)]
%%
if 1
x_iter = sol.W(ind_x,:);
% x_iter =x_iter(1:d+1:end,:);

x_iter1= x_iter(1:n_x:end,:);
x_iter2= x_iter(2:n_x:end,:);
figure
plot(x_iter1,x_iter2)

legend_str = {};
for ii = 1:size(x_iter2,2)
    legend_str  = [legend_str ; ['iter ' num2str(ii)]];
end
legend(legend_str);
end



