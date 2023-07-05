close all
clear all
import casadi.*

problem_options = NosnocProblemOptions();
problem_options.n_s = 2;
problem_options.dcs_mode = 'Stewart';
problem_options.print_level = 3;

model = NosnocModel();
x1 = SX.sym('x1');
model.x = x1;
model.c = [x1];
model.S = [-1; 1];

f_11 = 3;
f_12 = 1;
model.F = [f_11, f_12];
model.x0 = -1;

model.T = pi/4;
problem_options.N_stages = 1;
problem_options.N_finite_elements = 2;
model.dims.n_s = 2;

model.verify_and_backfill(problem_options);
model.generate_variables(problem_options);
model.generate_equations(problem_options);
problem_options.preprocess();

mpcc = NosnocMPCC(problem_options, model.dims, model);

solver_options = NosnocSolverOptions();
solver_options.preprocess();

%solver = NosnocSolver(model, settings);
%[results,stats] = solver.solve();
