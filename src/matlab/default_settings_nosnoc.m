function [default_settings] = default_settings_nosnoc()

% TODO: add descprition of function the list of default settings and their explanation.
%% General
solver_name = 'nosnoc_solver';
use_fesd = 1;
casadi_symbolic_mode = 'MX';
%% IRK anf FESD Settings
n_s = 2;                            % Number of IRK stages
irk_scheme = 'radau';     % RK scheme
irk_representation = 'integral'; % are the IRK equations in differential from (derivative at stages are uknowns in the equations) or in integral form (state values are unkwnowns at stage points)
lift_irk_differential = 1; % if differential mode is used, introduce new variables for intermediate stage values X_ki.
cross_comp_mode = 3;
gamma_h = 1;


pss_mode = 'Stewart'; % possible options: Stewart and Step


% initalization - Stewart
lp_initalization = 0;
initial_theta = 0;
initial_lambda = 0;
initial_mu = 0;
% initalization - Step
initial_alpha = 1;
initial_lambda_0 = 1;
initial_lambda_1 = 1;
initial_beta = 1;
initial_gamma = 1;

pss_lift_step_functions = 1; % lift the multilinear terms in the step functions;
n_depth_step_lifting = 2; % it is not recomended to change this (increase nonlinearity and harms convergenc)

couple_across_stages = 1;
list_of_all_rk_schemes = {'radau','legendre','Radau-IIA','Gauss-Legendre','Radau-I','Radau-IA',...
                           'Lobatto-III','Lobatto-IIIA','Lobatto-IIIB','Lobatto-IIIC',...
                           'Explicit-RK'};
%% General NLP/OCP Settings
g_ineq_constraint = 0; % is nonlinear path constraint present (by default evaluated only on control grid points)
g_ineq_at_fe = 0; % evaluate nonlinear path constraint at every finte element boundary
g_ineq_at_stg = 0; % evaluate nonlinear path constraint at every stage 

% x_box_at_fe = 1; % evaluate box constraint for diff states at every finite element boundary point
% x_box_at_stg = 1; % evulate box constraint for diff states at every stage point. (is set to zero per default in differential irk mode, as it becomes a linear instead of box constraint)

terminal_constraint = 0;
time_optimal_problem = 0;
simple_v0_guess = 0;
%% MPCC and Homotopy Settings	
comp_tol = 1e-14;
mpcc_mode = 5;
objective_scaling_direct = 1;
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
% sigma_0 = sigma_scale*rho_scale*exp(-rho_lambda*rho_0);
   
nonlinear_sigma_rho_constraint = 1;
convex_sigma_rho_constraint = 0;

ratio_for_homotopy_stop = 0.75;

polishing_step = 0; % huristic for fixing active set, yet exerimental, not recommended to use.
polishing_derivative_test = 0; % check in sliding mode also the derivative of switching functions

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
	
%% Time-Settting & Time-Freezing
time_freezing  = 0;
time_freezing_inelastic = 0;
time_optimal_problem = 0;
time_rescaling = 0;
use_speed_of_time_variables = 1;
local_speed_of_time_variable = 0;
stagewise_clock_constraint = 1;
s_sot0 = 1;
s_sot_max =	25;
s_sot_min =	1/s_sot_max;
impose_terminal_phyisical_time = 1;
T_final_max = 1e2;
T_final_min = 0;
time_freezing_reduced_model = 0; % analytic reduction of lifter formulation, less algebraic variables (experimental)
time_freezing_hysteresis = 0; % do not do automatic time freezing generation for hysteresis, it is not supported yet.
%% Verbose
print_level = 3;
opts_ipopt.ipopt.print_level = 0;
opts_ipopt.print_time = 0;
opts_ipopt.ipopt.sb= 'yes';

%% IPOPT Settings
opts_ipopt.verbose = false;
opts_ipopt.ipopt.max_iter = 500;
% opts_ipopt.ipopt.print_level = 5;
opts_ipopt.ipopt.bound_relax_factor = 0;
tol_ipopt = 1e-16;
opts_ipopt.ipopt.tol = tol_ipopt;
opts_ipopt.ipopt.dual_inf_tol = tol_ipopt;
opts_ipopt.ipopt.dual_inf_tol = tol_ipopt;
opts_ipopt.ipopt.compl_inf_tol = tol_ipopt;
opts_ipopt.ipopt.mu_strategy = 'adaptive';
opts_ipopt.ipopt.mu_oracle = 'quality-function';
opts_ipopt = opts_ipopt;
%% Relxation of terminal constraint
relax_terminal_constraint = 0; %  0  - hard constraint, 1 - ell_1 , 2  - ell_2 , 3 - ell_inf;
relax_terminal_constraint_from_above = 0; 
rho_terminal = 1e2;

%% Ingerator Specific 
use_previous_solution_as_initial_guess = 0;
simulation_problem  = 0;

%% Misc
missing_settings_already_filled_in = 1;
there_exist_free_x0 = 0;
clear_ipopt_verbose = 0;
output_stage_values = 0;
time_freezing_model_exists = 0;


%% Save data into struct
names = who;

for ii = 1:length(names)
    eval([ 'default_settings.' names{ii} '=' names{ii} ';'])
end

end

