clear all
close all
import casadi.*
import vdx.*

T = 2.0;
N_sim = 100;
t_step = T/N_sim;
%% Define (uncontrolled for now) projected system
x = SX.sym('x', 2);
data.x = x;
data.lbx = [-inf;-inf];
data.ubx = [inf;inf];
data.x0 = [0;0.5];
data.u = [];
data.lbu = [];
data.ubu = [];
data.u0 = [];
data.c = [x(2);-x(2) - (x(1)+0.5)^2 + 1];
data.f_x = [x(2); -x(1)];
data.f_q = 0;
data.f_q_T = 0;

data.T = t_step;
data.N_stages = 1;
data.N_fe = 3;
data.n_s = 2;
data.irk_scheme = 'radau';

prob = InclusionProblem(data, struct);

prob.generate_constraints();

prob.create_solver(struct);
