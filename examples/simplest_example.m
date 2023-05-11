close all
clear all
import casadi.*

settings = NosnocOptions();
settings.n_s = 2;
settings.dcs_mode = 'Step';
settings.mpcc_mode = MpccMode.Scholtes_eq;
settings.print_level = 3;

x1 = SX.sym('x1');
model.p_time_var = SX.sym('t_var',1);
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

solver = NosnocSolver(model, settings);
[results,stats] = solver.solve();
