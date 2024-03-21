clear all
close all
import casadi.*
import vdx.*

T = 0.35;
R = 3.5;
%% Define (uncontrolled for now) projected system
x1 = SX.sym('x1', 2);
x2 = SX.sym('x2', 2);
x = [x1;x2];
data.x = x;
data.lbx = [-inf;-inf;-inf;-inf];
data.ubx = [inf;inf;inf;inf];
data.x0 = [-25;-25;-15;-15];
u1 = SX.sym('u1', 2);
u2 = SX.sym('u2', 2);
data.u = [u1;u2];
data.lbu = [-100/sqrt(2);-100/sqrt(2);-60/sqrt(2);-60/sqrt(2)];
%data.ubu = [100/sqrt(2);100/sqrt(2);60/sqrt(2);60/sqrt(2)];
data.ubu = [100/sqrt(2);100/sqrt(2);0;0];
data.u0 = [0;0;0;0];
data.c = [norm_2(x2-x1)-2*R];
data.f_x = [u1;u2];
data.f_q = 0;%0.01*norm_2(data.u)^2;
data.f_q_T = 0.5*(norm_2(x)^2);

data.T = T;
data.N_stages = 10;
data.N_fe = 3;
data.n_s = 2;
data.irk_scheme = 'radau';

prob = InclusionProblem(data, struct);

prob.generate_constraints();

prob.create_solver(struct);


prob.w.x(0,0,0).init = data.x0;
prob.w.x(0,0,0).lb = data.x0;
prob.w.x(0,0,0).ub = data.x0;
homotopy(prob);
