clear all; 
close all; 
clc;
delete robot_ocp.mat
delete robot_ocp.gif
delete log_robot.txt

%% robot scene description
scenario.constant_inertia_matrix = 0;
scenario.relax_control_bounds = 0;
scenario.general_inequality_constraints = 1;
scenario.save_figure = 1;
scenario.filename = 'robot_ocp';
scenario.u0 = [0;0];
scenario.polishing_penalty_iteration = 0;

%% auxiliary dynamics and friction
scenario.a_n = 200;
scenario.mu = 0.8;
%% obstacles
scenario.hole_constraint = 1;
scenario.q_target = [3;0.4;0;0];
scenario.n_holes = 3;
scenario.xc_vec = [0.5 1.5 2.5];
scenario.zc_vec = [0 0 0];
scenario.height_vec = [0.1 0.1 0.1];
scenario.width_vec = [0.5 0.5 0.5];

%% Default settings NOSNOC
[settings] = NosnocOptions();
settings.print_level = 3;
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 2;
settings.pss_lift_step_functions = 0;
%% homotopy settings
settings.cross_comp_mode = 3;
settings.opts_ipopt.ipopt.max_iter = 2500;
settings.N_homotopy = 5;
settings.homotopy_update_rule = 'superlinear';
settings.opts_ipopt.ipopt.tol = 1e-6;
settings.opts_ipopt.ipopt.acceptable_tol = 1e-6;
settings.opts_ipopt.ipopt.acceptable_iter = 3;
settings.comp_tol = 1e-10;
settings.opts_ipopt.ipopt.linear_solver = 'ma57';

%% time-freezing
settings.s_sot_max = 5;
settings.s_sot_min = 1;
settings.rho_sot = 0.01;
settings.time_freezing = 1;
settings.pss_lift_step_functions = 1;
settings.stagewise_clock_constraint = 0;

%% Discretization
model.T = 2.5;
model.N_stages = 25;
model.N_FE = 3;

%% call function to discretize, solve and plot results
[results,stats,model] = solve_robot_ocp(model,settings,scenario);
