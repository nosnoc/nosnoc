% Example of a small MPCC example, solved via the nonsnoc MPCC solver mpccsol. 
% The origin is C stationary, (1,0) is S-stationary
close all;
clear;
import casadi.*;
import nosnoc.solver.*;

mpccsol_opts = nosnoc.solver.Options();  


x1 = SX.sym('x1');
x2 = SX.sym('x2');
% parameters
p = SX.sym('p');
x = [x1;x2];
f = (x1-1)^2+x2^3+x2^2+p;
x0 = [0;0];
p0 = 0;
lbx = [-inf;-inf];
ubx = [inf;inf];

mpcc = vdx.problems.Mpcc();

mpcc.w.x(1) = {x1};
mpcc.w.x(2) = {x2};
mpcc.p.p = {p, p0};
mpcc.G.x = {x1};
mpcc.H.x = {x2};
mpcc.f = f;

mpcc_method1 = 'scholtes_ineq';

mpccsol_opts.homotopy_steering_strategy = "ELL_1";
% the parameter relaxation/smoothing parameter s is minmized in the objective, its penality (1/sigma) is steered outside, cf. Table 1 in https://arxiv.org/pdf/2312.11022.pdf 
% mpccsol_opts.homotopy_steering_strategy = "DIRECT"; 
% the parameter relaxation/smoothing parameter sigma is steered outside the optimization , cf. Table 1 in https://arxiv.org/pdf/2312.11022.pdf

mpccsol_opts.calculate_stationarity_type = 1;

solver = mpccsol('generic_mpcc', mpcc_method1, mpcc, mpccsol_opts);


mpcc_results = solver('x0', x0,...
    'lbx', lbx,...
    'ubx', ubx,...
    'p',p0);

disp(mpcc_results.x)


%% Check stationarity of non optimal points
%w_optimal = mpcc_results.x;
%[w_polished, res_out, stat_type, n_biactive] = solver.check_multiplier_based_stationarity(w_optimal);
% Try calculating b-stationarity with polished result:
%[solution, improved_point, b_stat] = solver.check_b_stationarity(w_polished);

%fprintf('-------------------------------\n');
%w_non_optimal = [0;0]; % - C and not B
%[w_polished, res_out, stat_type, n_biactive] = solver.check_multiplier_based_stationarity(w_non_optimal);
% Try calculating b-stationarity with polished result:
%[solution, improved_point, b_stat] = solver.check_b_stationarity(w_polished);
