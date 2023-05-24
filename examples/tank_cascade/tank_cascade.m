function [model] = tank_cascade()
import casadi.*
model = NosnocModel();
%% Inital value
model.x0 = [0.1;0.1;0.1];
model.u0 = [0;0;0;0]; % guess for control variables
%% Variable defintion
% differential states
L1 = SX.sym('L1');
L2 = SX.sym('L2');
L3 = SX.sym('L3');
x = [L1;L2;L3];
model.x = x;
n_x = 3;

% lower and upper bounds
model.lbx = 1e-3*ones(n_x,1);
model.ubx = inf*ones(n_x,1);

%% Control
wu0 = SX.sym('wu0');
wu1 = SX.sym('wu1');
wu2 = SX.sym('wu2');
wu3 = SX.sym('wu3');
u = [wu0;wu1;wu2;wu3];
model.u = u;
n_u = length(model.u);
% Guess and Bounds
model.lbu  = 0.25*ones(n_u,1);
model.ubu  = 1.25*ones(n_u,1);
% Parameters
H1 = -0.5;
H2 = -0.5;
H3 = -0.5;
k0  = 0.1;
k1  = 0.1;
k2  = 0.1;
k3  = 0.1;
% Auxliary functions
A1 = 50;
A2 = 50; 
A3 = 50;
%% Switching Functions
% every constraint funcion corresponds to a simplex (note that the c_i might be vector valued)
c1 = -(L2-H2);
c2 = -(L3-H3);
% sign matrix for the modes
S1 = [1;-1];
S2 = [1;-1];
% c = [c1;c2];
model.c = {c1,c2};
model.S = {S1,S2};

F_input = 0;
f11 = [wu0*k0-wu1*k1*sqrt(L1-L2-H2);wu1*k1*sqrt(L1-L2-H2);-wu3*k3*sqrt(L3)];
f12 = [wu0*k0-wu1*k1*sqrt(L1);wu1*k1*sqrt(L1);-wu3*k3*sqrt(L3)];
% for c2, h2
f21 = [0;-wu2*k2*sqrt(L2-L3-H3);wu2*k2*sqrt(L2-L3-H3)];
f22 = [0;-wu2*k2*sqrt(L2);wu2*k2*sqrt(L2)];
% in matrix form
F1 = [f11 f12];
F2 = [f21 f22];
model.F = {F1,F2};

%% Objective
model.f_q = (x-ones(n_x,1)*0.75)'*(x-ones(n_x,1)*0.75);
model.f_q_T = 0;

end
