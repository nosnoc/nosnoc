
clear all
clc
close all
import casadi.*

%% Info
% This is an example from Stewart, D.E., 1996. A numerical method for friction problems with multiple contacts. The ANZIAM Journal, 37(3), pp.288-308.
% It considers 3 independent switching functions and it demonstrates the
% generalization of the FESD scheme to independet subystems (sum of Filippov systems)
%% discretization parameters
N_sim = 85;
T_sim = 12;
N_FE = 3; % (per integration step)

%% NOSNOC settings
% load general nosnoc problem options
problem_options = nosnoc.Options();
% load nosnoc mpcc solver options
solver_options = nosnoc.solver.Options();
solver_options.N_homotopy = 10;
solver_options.use_previous_solution_as_initial_guess = 1; 
% (this is more of an integrator options and should thus live in problem options?)

problem_options.n_s = 2;
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.rk_representation= 'differential';
problem_options.dcs_mode = 'Stewart';
problem_options.N_finite_elements = N_FE;
problem_options.N_sim = N_sim;
problem_options.T_sim = T_sim;
problem_options.cross_comp_mode = 3;

% create nosnoc pss model object
model = nosnoc.model.Pss();
model.x0 = [-1;1;-1;-1;1;1;0];
% differential states
q1 = SX.sym('q1');
q2 = SX.sym('q2');
q3 = SX.sym('q3');
v1 = SX.sym('v1');
v2 = SX.sym('v2');
v3 = SX.sym('v3');
t = SX.sym('t');
q = [q1;q2;q3];
v = [v1;v2;v3];
model.x = [q;v;t];
c1 = v1;
c2 = v2;
c3 = v3;
% sign matrix for the modes
S = [1;-1]; % same for all three subsystems
% discrimnant functions
model.S = {S,S,S};
model.c = {c1,c2,c3};
F_external = 0; % external force, e.g., control
F_input = 10; % variable force exicting
f_base = [v1;...
    v2;...
    v3;...
    (-q1)+(q2-q1)-v1;...
    (q1-q2)+(q3-q2)-v2;...
    (q2-q3)-v3+F_input*cos(pi*t);...
    ]/6;

f_base = [v1;...
    v2;...
    v3;...
    (-q1)+(q2-q1)-v1;...
    (q1-q2)+(q3-q2)-v2;...
    (q2-q3)-v3+F_external+F_input*(1*0+1*cos(pi*t));...
    1];
%
% for c1, h1,
f_11 = f_base+[0;0;0;-0.3;0;0;0];
f_12 = f_base+[0;0;0;+0.3;0;0;0];
% for c2, h2
f_21 = [0;0;0;0;-0.3;0;0];
f_22 = [0;0;0;0;0.3;0;0];
% for c3, h3
f_31 = [0;0;0;0;0;-0.3;0];
f_32 = [0;0;0;0;0;0.3;0];
% unfold_struct(model,'base');
% in matrix form
F1 = [f_11 f_12];
F2 = [f_21 f_22];
F3 = [f_31 f_32];

model.F = {F1 F2 F3};

integrator = nosnoc.Integrator(model, problem_options, solver_options);
[t_grid, x_res] = integrator.simulate();

integrator.dcs
integrator.discrete_time_problem.dcs
%%
figure
subplot(211)
plot(t_grid,x_res(1,:));
hold on
plot(t_grid,x_res(2,:));
plot(t_grid,x_res(3,:));
xlabel('$t$','interpreter','latex');
ylabel('$q(t)$','interpreter','latex');
legend({'$q_1(t)$','$q_2(t)$','$q_3(t)$'},'interpreter','latex');
grid on
subplot(212)
plot(t_grid,x_res(4,:));
hold on
plot(t_grid,x_res(5,:));
plot(t_grid,x_res(6,:));
xlabel('$t$','interpreter','latex');
ylabel('$v(t)$','interpreter','latex');
legend({'$v_1(t)$','$v_2(t)$','$v_3(t)$'},'interpreter','latex');
grid on
