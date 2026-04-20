%% Cardinality optimization

% solve random cardinality optimization problem, using the foromulation
% from:
% Complementarity Formulations of \ell_0-norm Optimization Problems
% Mingbin Feng, John E. Mitchell, Jong-Shi Pang, Xin Shen, Andreas W achter

clear all; clc; close all
import casadi.*
import nosnoc.solver.*;


%% Problem formulation
rng(1,"twister");
n = 100; % number of variables
m = 10;
A = rand(m,n);
b = rand(m,1);
c = 1-2*rand(n,1);
Q = rand(n,n);
e = ones(n,1);
gamma = 1e3;

x = SX.sym('x',n);
x_p = SX.sym('x_p',n);
x_n = SX.sym('x_n',n);
xi = SX.sym('xi',n);
w = [x;x_p;x_n;xi];


f = 0.5*x'*Q*x+c'*x+gamma*e'*(e-xi)+sin(x'*x);

% f = e'*(e-xi);

g = [A*x-b; x-(x_p-x_n)];
x0 = ones(4*n,1);
G = [xi;x_p];
H = [x_p+x_n;x_n];
lbx = -inf(4*n,1);
lbx(n+1:end) = 0;
ubx = 1e2*ones(4*n,1);
ubx(end-n:end) = 1;
lbg = [zeros(m,1);zeros(n,1)];
ubg = [inf(m,1);zeros(n,1)];
p = [];
p0 = [];

mpec = struct('w', w,'f',f,'g',g,'G',G,'H',H,'p',p);
solver_initalization = struct('x0', x0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg,'p',p0);
%%  mpecopt settings
solver_settings = mpecopt.Options();
solver_settings.settings_lpec.lpec_solver = "Gurobi";
% solver_settings.rho_TR_phase_i_init = 1e-1;
solver_settings.settings_lpec.lpec_solver = "Highs_casadi";

% solver_settings.rho_TR_phase_ii_init = 1e-1;
% solver_settings.relax_and_project_homotopy_parameter_steering = "Direct";
% solver_settings.initialization_strategy = "FeasibilityEllInfGeneral";
%solver = mpecopt.Solver(mpec, solver_settings);
%%
plugin = 'ccopt';
sqpec_opts = nosnoc.sqpec.Options();
sqpec_opts.max_iter = 20;
sqpec_opts.tol = 1e-6;
% sqpec_opts.tol = 1e-6;
sqpec_opts.qpec_solver_options = solver_settings;
sqpec_opts.qpec_solver_options.verbose_solver = false;
sqpec_opts.qpec_solver_options.verbose_summary = false;
% mpccsol_opts.qpec_solver_options.print_level = 0;
% mpccsol_opts.qpec_solver_options.homotopy_steering_strategy = "DIRECT";
% mpccsol_opts.qpec_solver_options.assume_lower_bounds = true;
mpccsol_opts.rho_d = 0;
sqpec_opts.discard_constraints_in_hessian = false;
% solver_initalization.x0 = result_mpecopt.x;

% mpccsol_opts = nosnoc.reg_homotopy.Options();
% plugin = 'reg_homotopy';
opts.casadi_opts = struct();
opts.casadi_opts.madmpec.bound_relax_factor = 0.0;
%opts.casadi_opts.madmpec.linear_solver = 'CUDSSSolver';
solver = mpccsol('generic_mpcc', plugin, mpec, opts);
mpcc_results = solver(solver_initalization);
w_opt_sqpec = mpcc_results.x;
f_opt_sqpec = mpcc_results.f;
%%

x_opt_sqpec = w_opt_sqpec(1:n);
cardinality_sqpec  = sum(heaviside(abs(x_opt_sqpec)-1e-3));
