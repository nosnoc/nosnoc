% Example 3 from Scheel, Holger, and Stefan Scholtes. "Mathematical programs with complementarity constraints: Stationarity, optimality, and sensitivity." Mathematics of Operations Research 25.1 (2000): 1-22.
% The solution is M or C stationary, and B stationary, but not S stationary. MPCC-LICQ is violated.
close all;
clear;
clc;
import casadi.*;
import nosnoc.solver.*;

% set options
mpccsol_opts = nosnoc.reg_homotopy.Options();  
mpccsol_opts.homotopy_steering_strategy = "Direct";
% mpccsol_opts.homotopy_steering_strategy = "ELL_1";
mpccsol_opts.homotopy_steering_strategy = "ELL_INF";


mpccsol_opts.complementarity_tol = 1e-10;

% define MPCC
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1;x2];
f = (x1-1)^2+(x2+1)^2;

x0 = [0;0];
lbx = [0;0];
ubx = [inf;inf];

g = [x2^2];
lbg = [-inf];
ubg = [0];

G = -x1;
H = x2;

mpcc.x = x;
mpcc.g = g;
mpcc.G = x1;
mpcc.H = x2;
mpcc.f = f;

solver = mpccsol('solver_mpecc', 'reg_homotopy', mpcc, mpccsol_opts);
mpcc_results = solver('x0', x0,'lbx', lbx,'ubx', ubx);
disp(mpcc_results.x)
%% Check stationarity points 

% Try calculating b-stationarity before getting polished result:
[solution, improved_point, b_stat] = solver.check_b_stationarity(mpcc_results.x);
% get polished solution and multiplier based stationarity type with a lifted TNLP.
[w_polished, res_out, stat_type, n_biactive] = solver.check_multiplier_based_stationarity(mpcc_results.x);
% Try calculating b-stationarity with polished result:
[solution, improved_point, b_stat] = solver.check_b_stationarity(w_polished);

[w_polished, res_out, stat_type, n_biactive] = solver.check_multiplier_based_stationarity([0;0]);
