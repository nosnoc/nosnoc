% Example with general nonlinear constraints
close all;
clear;
import casadi.*;
import nosnoc.solver.mpccsol;

mpccsol_opts = nosnoc.solver.Options();
mpccsol_opts.opts_casadi_nlp.ipopt.mu_strategy = 'adaptive';
mpccsol_opts.opts_casadi_nlp.ipopt.mu_oracle = 'quality-function';

x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
x4 = SX.sym('x4');
x5 = SX.sym('x5');
x6 = SX.sym('x6');
x7 = SX.sym('x7');
x8 = SX.sym('x8');
w = [x1;x2;x3;x4;x5;x6;x7;x8];
p = SX.sym('p');
f = (x1-5)^2+(2*x2+1)^2+x4^4+x5^3+p^2;
g = [2*(x2-1)-1.5*x1+x3-0.5*x4+x5
    x1^2+x2^2+x3^2-10
    -x1+0.5*x2+4-x7^2;...
    -x1-x2^2+1-x8^3;...
    x1^2+x2^2+x3^2+x4^2-16;...
    ];
G = [x6;x7;x8];
H = [x3;x4;x5];

w0 = 5*ones(8,1);
lbw = [-2;-2;zeros(6,1)];
ubw = [1;2;3;inf*ones(5,1)];
lbg = [zeros(4,1);-inf;];
ubg = zeros(5,1);

p0 = 0;

mpcc_struct.x = w;
mpcc_struct.g = g;
mpcc_struct.p = p;
mpcc_struct.G = G;
mpcc_struct.H = H;
mpcc_struct.f = f;

solver = mpccsol('generic_mpcc', 'scholtes_ineq', mpcc_struct, mpccsol_opts);

mpcc_results = solver('x0', w0,...
    'lbx', lbw,...
    'ubx', ubw);

disp(mpcc_results.x)
disp(mpcc_results.f)
