function [results,stats] = test_generic_mpcc(mpcc_method, homotopy_steering_strategy)
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
lbx = [0;0];
ubx = [inf;inf];

mpcc_struct.x = x;
mpcc_struct.g = [];
mpcc_struct.p = p;
mpcc_struct.G = x1;
mpcc_struct.H = x2;
mpcc_struct.f = f;

mpcc_method = mpcc_method;

mpccsol_opts.homotopy_steering_strategy = homotopy_steering_strategy;
mpccsol_opts.relaxation_strategy = mpcc_method;


solver = mpccsol('generic_mpcc', 'reg_homotopy', mpcc_struct, mpccsol_opts);

results = solver('x0', x0,...
    'lbx', lbx,...
    'ubx', ubx,...
    'p',p0);

stats = solver.stats;
end
