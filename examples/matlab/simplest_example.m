close all
clear all
import casadi.*
[settings] = default_settings_nosnoc();
settings.n_s = 2;
settings.pss_mode = 'Stewart';
settings.mpcc_mode = 'Scholtes_ineq';
settings.step_equilibration = StepEquilibrationMode.heuristic_diff;
settings.initial_lambda = 1;
settings.initial_theta = 1;
settings.initial_mu = 1;
settings.print_level = 3;
x1 = SX.sym('x1');
model.x = x1;
model.c = [x1];
model.S = [-1; 1];

f_11 = 3;
f_12 = 1;
model.F = [f_11, f_12];
model.x0 = -1;

model.T = pi/4;
model.N_stages = 1;
model.N_finite_elements = 2;

[results,stats,model,settings] = nosnoc_solver(model,settings);
