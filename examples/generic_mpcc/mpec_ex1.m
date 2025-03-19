% Example of a small MPCC example, solved via the nonsnoc MPCC solver mpccsol. 
% The origin is C stationary, (1,0) is S-stationary
close all;
clear;
import casadi.*;
import nosnoc.solver.*;

%%
mpccsol_opts = nosnoc.reg_homotopy.Options();  % Options for the relxation homotopy solvers in nosnoc
mpecopt_opts = mpecopt.Options(); % Options for the external active set mpecopt solver (see https://github.com/nosnoc/mpecopt)

%% Define MPCC problem
x1 = SX.sym('x1');
x2 = SX.sym('x2');
% parameters
p = SX.sym('p');
x = [x1;x2];
f = (x1-1)^2+x2^3+x2^2+p;
x0 = [0;0];
p0 = 0;
lbx = [0;0];
ubx = [inf;inf];

mpcc_data.x = x;
mpcc_data.g = x1+x2;
mpcc_data.p = p;
mpcc_data.G = x1;
mpcc_data.H = x2;
mpcc_data.f = f;

solver_initialization = struct('x0', x0,'lbx', lbx,'ubx', ubx,'p',p0,'lbg',-inf,'ubg',10);

mpcc_method1 = 'scholtes_ineq';
mpcc_method2 = 'kanzow_schwartz_ineq';
mpcc_method3 =  nosnoc.reg_homotopy.MpccMethod.STEFFENSEN_ULBRICH_INEQ;

mpccsol_opts.homotopy_steering_strategy = "ELL_1";
% the parameter relaxation/smoothing parameter s is minmized in the objective, its penality (1/sigma) is steered outside, cf. Table 1 in https://publications.syscop.de/Nurkanovic2024b.pdf
% mpccsol_opts.homotopy_steering_strategy = "DIRECT"; 
% the parameter relaxation/smoothing parameter sigma is steered outside the optimization , cf. Table 1 in https://publications.syscop.de/Nurkanovic2024b.pdf

%% Create several different solvers
mpccsol_opts.calculate_stationarity_type  = true;
mpccsol_opts.relaxation_strategy = mpcc_method1;
solver = mpccsol('scholtes_reg', 'reg_homotopy', mpcc_data, mpccsol_opts);

mpccsol_opts.relaxation_strategy = mpcc_method2;
solver_ks = mpccsol('kanzow_schwartz_solver', 'reg_homotopy', mpcc_data, mpccsol_opts);

mpccsol_opts.relaxation_strategy = mpcc_method3;
solver_su = mpccsol('su_solver', 'reg_homotopy', mpcc_data, mpccsol_opts);

% Mpec opt solver object
solver_mpecopt = mpccsol('generic_mpcc', 'mpecopt', mpcc_data, mpecopt_opts); % call mpecopt via nosnoc
% solver_mpecopt = mpecopt.Solver(mpcc_data,mpecopt_opts); % call mpecopt directly via its native interface

%% Test different solvers
mpcc_results = solver(solver_initialization);
disp(mpcc_results.x)

mpcc_results = solver_ks('x0', x0,'lbx', lbx,'ubx', ubx,'p',p0,'lbg',-inf,'ubg',10); % solver can be called by passing a struct (as above) or directly like CasADi's nlpsol
disp(mpcc_results.x)

mpcc_results = solver_su(solver_initialization);
solver_ks.compute_constraint_violation()
solver_su.stats
disp(mpcc_results.x)

mpcc_results = solver_mpecopt(solver_initialization);
disp(mpcc_results.x)

%% Check stationarity of optimal and non-optimal points
w_optimal = mpcc_results.x;
[w_polished, res_out, stat_type, n_biactive] = solver.check_multiplier_based_stationarity(w_optimal);
% Try calculating b-stationarity with polished result:
[solution, improved_point, b_stat] = solver.check_b_stationarity(w_polished);

fprintf('-------------------------------\n');
w_non_optimal = [0;0]; % - C and not B
[w_polished, res_out, stat_type, n_biactive] = solver.check_multiplier_based_stationarity(w_non_optimal);
% Try calculating b-stationarity with polished result:
[solution, improved_point, b_stat] = solver.check_b_stationarity(w_polished);
