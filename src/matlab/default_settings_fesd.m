function [default_settings] = default_settings_fesd()

% TODO: add descprition of function the list of default settings and their explanation.
%% General
solver_name = 'solver_fesd';

%% IRK anf FESD Settings
d = 2;                            % Degree of interpolating polynomial
collocation_scheme = 'radau';     % Collocation scheme: radau or legendre
N_stages = 10;                     %   
N_finite_elements = 3;            % Number of finite elments on a stage/control interval 

use_fesd = 1;
T = 2;
h = T/N_stages;
h_k = h/N_finite_elements;
fesd_complementartiy_mode = 1; % is replaced by cross_complementarity_mode 
cross_complementarity_mode = 3;
gamma_h = 1;
% initalization
lp_initalization = 0;
initial_theta = 0;
initial_lambda = 0;
initial_mu = 0;

% an experimental variable. should be always one if no equdistant control
% grid is
couple_across_stages = 1;
%% General NLP/OCP Settings
general_nonlinear_constraint = 0;
general_nonlinear_constraint_at_collocation_points = 0;
terminal_constraint = 0;
time_optimal_problem = 0;

%% MPCC and Homotopy Settings	
comp_tol = 1e-14;
mpcc_mode = 5;
objective_scaling_direct = 1;
% pointwise_or_integral = 1;
sigma_0	= 1;
sigma_N	= comp_tol;
kappa =	0.1;
N_homotopy = ceil(abs(log(sigma_N/sigma_0)/log(kappa)));
s_elastic_0	= 1;
s_elastic_max =	1e1;
s_elastic_min = 0;
store_all_homotopy_iterates = 1;

% Default settings for the barrier tuned penalty/slack variables for mpcc modes 8 do 10.
rho_penalty = 1e1;
sigma_penalty = 0;

rho_lambda = 1;
rho_scale = 30;

sigma_scale = 0.1;

rho_min = 0.1;
rho_max = (log(rho_scale)-log(1e-16))/rho_lambda;
  
rho_0 = max(rho_min,0.5);
sigma_0 = sigma_scale*rho_scale*exp(-rho_lambda*rho_0);
   
nonlinear_sigma_rho_constraint = 1;
convex_sigma_rho_constraint = 0;

%% Step equilibration	
regularize_h =1;
rho_h = 1;
delta_h_regularization = 0;
piecewise_equidistant_grid	= 0;
piecewise_equidistant_grid_sigma =	1;
piecewise_equidistant_grid_slack_mode = 0;

step_equilibration = 0;
step_equilibration_mode = 1;
step_equilibration_penalty = 0.1;  %(rho_h in step_equilibration modde 1, as qudratic penalty)
treat_step_equilibration_via_mpcc = 0;
step_equilibration_sigma = 0.1;
heuristic_step_equilibration = 1;
heuristic_step_equilibration_mode = 1; % 1 standard (h_k - h), 2  - penalize delta h
% heuristic_step_equilibration_penalty = 0.1;  %(rho_h)
%% Multiple shooting type discretization	
equidistant_control_grid = 1;
	
%% Time-Settting	
time_freezing  = 0;
time_rescaling = 0;
use_speed_of_time_variables = 1;
local_speed_of_time_variable = 1;
stagewise_clock_constraint = 1;
s_sot0 = 1;
s_sot_max =	25;
s_sot_min =	1/s_sot_max;
impose_terminal_phyisical_time = 1;

%% IPOPT Settings
opts_ipopt.verbose = false;
opts_ipopt.ipopt.max_iter = 500;
opts_ipopt.ipopt.print_level = 5;
opts_ipopt.ipopt.bound_relax_factor = 0;
tol_ipopt = 1e-16;
opts_ipopt.ipopt.tol = tol_ipopt;
opts_ipopt.ipopt.dual_inf_tol = tol_ipopt;
opts_ipopt.ipopt.dual_inf_tol = tol_ipopt;
opts_ipopt.ipopt.compl_inf_tol = tol_ipopt;
opts_ipopt.ipopt.mu_strategy = 'adaptive';
opts_ipopt.ipopt.mu_oracle = 'quality-function';
opts_ipopt = opts_ipopt;

%% Ingerator Specific 
use_previous_solution_as_initial_guess = 0;

%% Misc
explicit_z00 = 0; %% TODO: remove this when clear, it is an option during developemtn
missing_settings_already_filled_in = 1;
there_exist_free_x0 = 0;
%% Save data into struct
names = who;

for ii = 1:length(names)
    eval([ 'default_settings.' names{ii} '=' names{ii} ';'])
end

end

