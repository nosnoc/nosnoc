% Example 3 from Scheel, Holger, and Stefan Scholtes. "Mathematical programs with complementarity constraints: Stationarity, optimality, and sensitivity." Mathematics of Operations Research 25.1 (2000): 1-22.
% The solution is M or C stationary, and B stationary, but not S stationary. MPCC-LICQ is violated.
close all;
clear;
clc;
import casadi.*;
import nosnoc.solver.*;

% set options
mpccsol_opts = nosnoc.solver.Options();  
% mpccsol_opts.homotopy_steering_strategy = "Direct";
% mpccsol_opts.homotopy_steering_strategy = "ELL_1";
mpccsol_opts.homotopy_steering_strategy = "ELL_INF"; 
mpccsol_opts.calculate_stationarity_type = 1;

mpccsol_opts.complementarity_tol = 1e-12;

% mpcc_method =  nosnoc.solver.MpccMethod.KADRANI;
mpcc_method =  nosnoc.solver.MpccMethod.SCHOLTES_INEQ;


% define MPCC
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
x = [x1;x2;x3];
f = x1+x2-x3;

x0 = [0;0;0];
lbx = [0;0;-inf];
ubx = [inf;inf;inf];

g = [4*x1-x3;...
     4*x2-x3];
lbg = [0;0];
ubg = [inf;inf];

G = x1;
H = x2;


mpcc.x = x;
mpcc.g = g;
mpcc.G = x1;
mpcc.H = x2;
mpcc.f = f;

solver = mpccsol('solver_mpecc', mpcc_method, mpcc, mpccsol_opts);

mpcc_results = solver('x0', x0,'lbx', lbx,'ubx', ubx);
disp(mpcc_results.x)
