close all
clear all
import casadi.*

settings = NosnocOptions();
settings.n_s = 2;
settings.dcs_mode = 'Stewart';
settings.mpcc_mode = MpccMode.ell_1_penalty;
settings.print_level = 3;

model = NosnocModel();
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
settings.N_stages = 1;
settings.N_finite_elements = 2;
model.dims.n_s = 2;



solver = NosnocSolver(model, settings);
[results,stats] = solver.solve();
