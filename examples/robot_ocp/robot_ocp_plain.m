clear all; close all; clc;
%% robot scene description
scenario.constant_inertia_matrix = 0;
scenario.relax_control_bounds = 0;
scenario.general_inequality_constraints = 1;
scenario.save_figure = 0;
scenario.filename = 'robot_ocp_plain';
scenario.u0 = [0;0];
scenario.polishing_penalty_iteration = 0;
%% auxiliary dynamics and friction
scenario.a_n = 200;
scenario.mu = 0.8;
%% obstacles  - no obastacles, just plain gaol reaching
scenario.hole_constraint = 0;
scenario.q_target = [3;0.4;0;0];
scenario.n_holes = 0;
scenario.xc_vec = [];
scenario.zc_vec = [];
scenario.height_vec = [];
scenario.width_vec = [];

%% Default settings NOSNOC
[settings] = NosnocOptions();
model = NosnocModel();
settings.use_fesd = 1;
settings.print_level = 4;
settings.irk_scheme = IRKSchemes.RADAU_IIA;
%% homotopy settings
settings.cross_comp_mode = 3;
settings.opts_casadi_nlp.ipopt.max_iter = 1500;
settings.N_homotopy = 10;
settings.opts_casadi_nlp.ipopt.tol = 1e-12;
settings.comp_tol = 1e-10;
settings.sigma_0 = 1;
settings.opts_casadi_nlp.ipopt.linear_solver = 'ma57';
settings.mpcc_mode = MpccMode.Scholtes_ineq;
%% time-freezing
settings.s_sot_max = 5;
settings.s_sot_min = 1;
settings.rho_sot = 0.01;
settings.time_freezing = 1;
settings.pss_lift_step_functions = 1;
settings.stagewise_clock_constraint = 0;

%% Discretization
problem_options.T = 2.5;
settings.N_stages = 20;
settings.n_s = 2;
settings.N_finite_elements = 3;
%% call function to discretize, solve and plot results
[results,stats,model] = solve_robot_ocp(model,settings,scenario);
