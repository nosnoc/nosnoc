% Example from Kirches, C., Larson, J., Leyffer, S., & Manns, P. (2022).
% Sequential linearization method for bound-constrained mathematical programs
% with complementarity constraints.
% SIAM Journal on Optimization, 32(1), 75-99. Section 5.3
close all;
clear;
import casadi.*;
import nosnoc.solver.mpccsol;

mpccsol_opts1 = nosnoc.solver.Options();
mpccsol_opts1.lift_complementarities = true;
mpccsol_opts1.calculate_stationarity_type = true;
mpccsol_opts2 = nosnoc.solver.Options();
mpccsol_opts2.lift_complementarities = true;
mpccsol_opts2.calculate_stationarity_type = true;
x = SX.sym('x',8);
x_0 = x(1:4);
x_1 = x(5:6);
x_2 = x(7:8);

rho = 2;
lambda1= 3.9375;
lambda2 = -6.5;
lambda3 = -0.25;
lambda4 = 2.5;

f = 0.5*((x_0(1)-x_0(3))^2+(x_0(2)-x_0(4))^2)...
    +lambda1*(-34+2*x_0(3)+8/3*x_0(4)+x_2(1))...
    -lambda2*(-24.25+1.25*x_0(3)+2*x_0(4)+x_2(2))...
    -lambda3*(x_1(1)+x_0(2)+x_0(3)-15)...
    +lambda4*(x_1(2)+x_0(1)-x_0(4)-15)...
    +0.5*rho*( (-34+2*x_0(3)+8/3*x_0(4)+x_2(1))^2 + (-24.25+1.25*x_0(3)+2*x_0(4)+x_2(2))^2+ (x_1(1)+x_0(2)+x_0(3)-15)^2+ (x_1(2)+x_0(1)-x_0(4)-15)^2 );


x0 = zeros(8,1);
% x0 = [  11.6667
%     3.2083
%    11.6667
%     3.2083
%     0.0000
%     5.2917
%     0.1424
%    -0.0000];

lbx = [-inf*ones(4,1);-inf*ones(4,1)];
ubx = inf*ones(8,1);

mpcc_struct.x = x;
mpcc_struct.g = [];
mpcc_struct.p = [];
mpcc_struct.G = x_1;
mpcc_struct.H = 2*x_2;
mpcc_struct.f = f;

% create two solvers with different mpcc relaxation methods: steffensen_ulbrich_eq and steffensen_ulbrich_ineq
mpccsol_opts1.relaxation_strategy = 'steffensen_ulbrich_eq';
solver_eq = mpccsol('generic_mpcc', 'reg_homotopy', mpcc_struct, mpccsol_opts1);
mpccsol_opts2.relaxation_strategy = 'steffensen_ulbrich_ineq';
solver_ineq = mpccsol('generic_mpcc', 'reg_homotopy', mpcc_struct, mpccsol_opts2);

mpcc_results = solver_eq('x0', x0,...
    'lbx', lbx,...
    'ubx', ubx);

disp(mpcc_results.x)
disp(mpcc_results.f)

mpcc_results = solver_ineq('x0', x0,...
    'lbx', lbx,...
    'ubx', ubx);

disp(mpcc_results.x)
disp(mpcc_results.f)
