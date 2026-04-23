% Example of a small MPCC example, solved via the nonsnoc MPCC solver mpccsol with sqpec plug in.
close all;
clear;
import casadi.*;
import nosnoc.solver.*;

use_vdx = 0;


plugin = 'sqpec';
%plugin = 'ccopt';
% plugin = 'reg_homotopy';
% plugin = 'mpecopt';

switch plugin
  case 'sqpec'
    mpccsol_opts = nosnoc.sqpec.Options();
    mpccsol_opts.max_iter = 5;
    mpccsol_opts.tol = 1e-6;
    mpccsol_opts.qpec_solver_options.print_level = 0;
  case 'reg_homotopy'
    mpccsol_opts = nosnoc.reg_homotopy.Options();
  case 'mpecopt'
    mpccsol_opts = mpecopt.Options();
  case 'ccopt'
    mpccsol_opts = struct();
    mpccsol_opts.casadi_opts = struct();
    mpccsol_opts.casadi_opts.madmpec = struct();
    mpccsol_opts.casadi_opts.madmpec.bound_relax_factor = 0.0;
    mpccsol_opts.casadi_opts.madmpec.print_level = 'TRACE';
    class(mpccsol_opts.casadi_opts.madmpec.bound_relax_factor)
end

x1 = SX.sym('x1');
x2 = SX.sym('x2');
% parameters
p = SX.sym('p');
x = [x1;x2];
f = (x1-1)^2+x2^3+x2^2;
g = x1.^4+x2.^4-p^2;
x0 = [0;2]; % needs to look around corner
            % x0 = [4;8];
p0 = 2;
% lbx = -inf(2,1);  % Todo; test if bounds are found
lbx = zeros(2,1);  % Todo; test if bounds are found
ubx = inf(2,1);
ubg = 0;
lbg = -inf;

if use_vdx
    mpec = vdx.problems.Mpcc();
    mpec.w.x(1) = {x1};
    mpec.w.x(2) = {x2};
    mpec.p.p = {p, p0};
    mpec.G.x = {x1};
    mpec.H.x = {x2};
    mpec.f = f;
    mpec.g.sym = g;
else
    mpec = struct();
    mpec.x = x;
    mpec.G = x1;
    mpec.H = x2 + x1;
    mpec.p = p;
    mpec.f = f;
    mpec.g = x1.^2+x2.^2-2^2;
end

%% Create solver object
solver = mpccsol('generic_mpcc', plugin, mpec, mpccsol_opts);
disp(solver)
% solve problem
solver_initalization = struct('x0', x0,'lbx', lbx,'ubx', ubx,'p',p0,'lbg',lbg,'ubg',ubg);
% mpcc_results = solver('x0', x0,'lbx', lbx,'ubx', ubx,'p',p0,'lbg',lbg,'ubg',ubg);
mpcc_results = solver(solver_initalization);
disp(mpcc_results.x)

