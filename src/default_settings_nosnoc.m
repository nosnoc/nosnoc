function [default_settings] = default_settings_nosnoc()

% TODO: add descprition of function the list of default settings and their explanation.
%% General
solver_name = 'nosnoc_solver';
use_fesd = 1;
casadi_symbolic_mode = 'SX';
%% Interface settings
real_time_plot = 0; % Plot solution trajectories during execution
%% IRK and FESD Settings
n_s = 2;                            % Number of IRK stages
irk_scheme = 'radau';     % RK scheme
irk_representation = 'integral'; % are the IRK equations in differential from (derivative at stages are uknowns in the equations) or in integral form (state values are unkwnowns at stage points)
lift_irk_differential = 1; % if differential mode is used, introduce new variables for intermediate stage values X_ki.
cross_comp_mode = 3;
gamma_h = 1;

pss_mode = 'Stewart'; % possible options: Stewart and Step

% initialization - Stewart
lp_initialization = 0;
initial_theta = 1;
initial_lambda = 1;
initial_mu = 1;
% initialization - Step
initial_alpha = 1;
initial_lambda_0 = 1;
initial_lambda_1 = 1;
initial_beta = 1;
initial_gamma = 1;

pss_lift_step_functions = 1; % lift the multilinear terms in the step functions;
n_depth_step_lifting = 2; % it is not recomended to change this (increase nonlinearity and harms convergenc), depth is number of multilinar terms to wich a lifting variables is equated to.

list_of_all_rk_schemes = {'radau','legendre','Radau-IIA','Gauss-Legendre','Radau-I','Radau-IA',...
                           'Lobatto-III','Lobatto-IIIA','Lobatto-IIIB','Lobatto-IIIC',...
                           'Explicit-RK'};
%% General NLP/OCP Settings
g_path_constraint = 0; % is nonlinear path constraint present (by default evaluated only on control grid points)
g_comp_path_constraint = 0;
g_path_at_fe = 0; % evaluate nonlinear path constraint at every finte element boundary
g_path_at_stg = 0; % evaluate nonlinear path constraint at every stage 

x_box_at_fe = 1; % evaluate box constraint for diff states at every finite element boundary point
x_box_at_stg = 1; % evaluate box constraint for diff states at every stage point. (is set to zero per default in differential irk mode, as it becomes a linear instead of box constraint)

terminal_constraint = 0;
time_optimal_problem = 0;
simple_v0_guess = 0;
%% MPCC and Homotopy Settings	
comp_tol = 1e-9;
mpcc_mode = 'Scholtes_ineq'; % 'direct', 'Scholtes_eq', 'Scholtes_ineq', 'ell_1_penalty', 'elastic_ineq', 'elastic_eq' , 'elastic_two_sided', 
% 'elastic_ell_1_ineq', 'elastic_ell_1_eq', 'elastic_ell_1_two_sided'
objective_scaling_direct = 1;
sigma_0	= 1;
sigma_N	= comp_tol;
homotopy_update_rule = 'linear'; % 'linear' sigma_k = homotopy_update_slope*sigma_N
    % 'superlinear' - sigma_k = max(sigma_N,min(homotopy_update_slope*sigma_k,sigma_k^homotopy_update_exponent));
homotopy_update_slope =	0.1;
homotopy_update_exponent = 1.5; % the exponent in the superlinear rule
N_homotopy = ceil(abs(log(sigma_N/sigma_0)/log(homotopy_update_slope)));
s_elastic_0	= 1;
s_elastic_max =	1e1;
s_elastic_min = 0;

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

%% Homotopy preprocess and polishing steps
h_fixed_iterations = 0;
h_fixed_max_iter = 1; % number of iterations that are done with fixed h in the homotopy loop
h_fixed_change_sigma = 1; % if this is on, do not update sigma and just solve on nlp with fixed h.
polishing_step = 0; % heuristic for fixing active set, yet exerimental, not recommended to use.
polishing_derivative_test = 0; % check in sliding mode also the derivative of switching functions
h_fixed_to_free_homotopy = 0; % start with large penaly for equidistant grid, end with variable equilibrated grid. 


%% Step equilibration
rho_h = 1;
step_equilibration = 'heuristic_mean'; % heuristic_mean, l2_relaxed, l2_relaxed_scaled, direct, direct_homotopy, direct_homotopy_lift
step_equilibration_sigma = 0.1; % slope at zero in rescaling the indicator function, nu_ki_rescaled = tanh(nu_ki/step_equilibration_sigma);

%% Multiple shooting type discretization	
equidistant_control_grid = 1;
	
%% Time-Setting & Time-Freezing
time_freezing  = 0;
time_freezing_inelastic = 0;
time_optimal_problem = 0;
time_rescaling = 0;
% for time optimal problems and equidistant control grids in physical time
use_speed_of_time_variables = 1;
local_speed_of_time_variable = 0;
stagewise_clock_constraint = 1;
impose_terminal_phyisical_time = 1;
s_sot0 = 1;
s_sot_max =	25;
s_sot_min =	1;
S_sot_nominal = 1;
rho_sot = 0;

T_final_max = 1e2;
T_final_min = 0;
time_freezing_reduced_model = 0; % analytic reduction of lifter formulation, less algebraic variables (experimental)
time_freezing_hysteresis = 0; % do not do automatic time freezing generation for hysteresis, it is not supported yet.
time_freezing_nonlinear_friction_cone = 1; % 1 - use nonlienar friction cone, 0 - use polyhedral l_inf approximation.

time_freezing_quadrature_state = 0; % make a nonsmooth quadrature state to integrate only if physical time is running
time_freezing_lift_forces = 0; % replace \dot{v} = M(q)^{-1}f(q,v,u) by dot{v} = z,  M(q)z - f(q,v,u) = 0; 

friction_exists = 0;
% exerimentla expert otpions
nonsmooth_switching_fun = 0; % experimental: use c = max(c1,c2) insetad of c = [c1;c2]  

%% Verbose
print_level = 3;

%% IPOPT Settings
opts_ipopt.ipopt.print_level = 0;
opts_ipopt.print_time = 0;
opts_ipopt.ipopt.sb = 'yes';
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

%% Relaxation of terminal constraint
relax_terminal_constraint = 0; %  0  - hard constraint, 1 - ell_1 , 2  - ell_2 , 3 - ell_inf;
relax_terminal_constraint_from_above = 0; 
rho_terminal = 1e2;
relax_terminal_constraint_homotopy = 0; % terminal penalty is governed by homotopy parameter

%% Integrator Specific 
use_previous_solution_as_initial_guess = 0;
simulation_problem  = 0;

%% Misc
there_exist_free_x0 = 0;
% missing_settings_already_filled_in = 1;
% clear_ipopt_verbose = 0;
% output_stage_values = 0;
time_freezing_model_exists = 0;


%% All NLP parameters
T_val = 1;
p_val = [sigma_0,rho_sot,rho_h,rho_terminal,T_val];

%% Save data into struct
names = who;

for ii = 1:length(names)
    eval([ 'default_settings.' names{ii} '=' names{ii} ';'])
end

end

