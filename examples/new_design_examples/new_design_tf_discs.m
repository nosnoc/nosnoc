% BSD 2-Clause License

% Copyright (c) 2022, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% This file is part of NOSNOC.

%% Manipulation of two discs
clear; close all; clc;
import casadi.*
import nosnoc.*
%%
filename = 'discs_manipulation.gif';
%% discretizatioon
N_stg = 15; % control intervals
N_FE = 2;  % integration steps per control interval
T = 2;
%% init
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
model = nosnoc.model.Cls();
%% settings
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 2;  % number of stages in IRK methods (TODO@Anton, other variations, e.g. Stewart, with n_s>1 converge nicely, this setting very slow)
problem_options.cross_comp_mode = "STAGE_STAGE";
problem_options.T = T;
problem_options.N_stages = N_stg;
problem_options.N_finite_elements = N_FE;

problem_options.time_freezing = 1;
problem_options.a_n = 10;
problem_options.relax_terminal_physical_time = ConstraintRelaxationMode.ELL_1;
problem_options.rho_terminal_physical_time = 1e5;
problem_options.use_numerical_clock_state = false;
problem_options.time_freezing_quadrature_state = true;

problem_options.time_freezing_Heaviside_lifting = true; 
problem_options.dcs_mode = "Stewart";

%% Solver settings
%solver_options.homotopy_update_rule = 'superlinear';
solver_options.homotopy_update_slope = 0.1;
solver_options.sigma_0 = 0.1;
solver_options.N_homotopy = 100;
solver_options.complementarity_tol = 1e-8;
solver_options.opts_casadi_nlp.ipopt.max_iter = 5e3;
solver_options.print_level = 3;
solver_options.warm_start_duals = true;
% IF HLS solvers for Ipopt installed (check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions) use the settings below for better perfmonace:
% solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';

%% model parameters
m1 = 2;
m2 = 1;
r1 = 0.3;
r2 = 0.2;

q10 = [-1; -2];
q20 = [-1;-1];
v10 = [0;0];
v20 = [0;0];

q_target1 = [-1; 1];
q_target2 = [0; 0];

x0 = [q10;q20;v10;v20];
ubx = [10; 10;10; 10; 5; 5; 5; 5];
lbx = -ubx;
u_max = [20;20];
u_min = -u_max;

u_ref = [0;0];
x_ref = [q_target1;q_target2;zeros(4,1)];
Q = diag([5;5;10;10;0*ones(4,1)]);
R = diag([0.1 0.1]);
Q_terminal = 100*Q;

%% Symbolic variables and bounds
q = SX.sym('q',4);
v = SX.sym('v',4);
u = SX.sym('u',2);

q1 = q(1:2);
q2 = q(3:4);

x = [q;v];
model.x = x;
model.u = u;
model.e = 0;
model.mu = 0.0;
model.x0 = x0;

model.M = diag([m1;m1;m2;m2]); % inertia/mass matrix;
model.f_v = [u;zeros(2,1)];

% gap functions
model.f_c = norm(q1-q2)^2-(r1+r2)^2;
model.dims.n_dim_contact = 2;

% box constraints on controls and states
model.lbu = u_min;
model.ubu = u_max;
model.lbx = lbx;
model.ubx = ubx;

%% Stage cost
model.f_q = (x-x_ref)'*Q*(x-x_ref)+ u'*R*u;
model.f_q_T = (x-x_ref)'*Q_terminal*(x-x_ref);

%% Call nosnoc solver
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();

%% read and plot results
x_res = ocp_solver.get('x');
h_res = ocp_solver.get('h');
p_res = x_res(1:4,:);
v_res = x_res(5:end,:);
p1 = x_res(1,:);
p2 = x_res(2,:);
p3 = x_res(3,:);
p4 = x_res(4,:);
v1 = x_res(5,:);
v2 = x_res(6,:);
v3 = x_res(7,:);
v4 = x_res(8,:);
t_opt = x_res(end,:);

frozen = logical([0, diff(t_opt) < 1e-3,]);
p_unfrozen = p_res;
v_unfrozen = v_res;
p_unfrozen(:, frozen) = [];
v_unfrozen(:, frozen) = [];
%% animation
t = linspace(0,2*pi,360).';t(end) = [];
pgon1 = polyshape(r1*cos(t), r1*sin(t));
pgon2 = polyshape(r2*cos(t), r2*sin(t));
facecolor1 = [0 0.4470 0.7410];
linecolor1 = facecolor1*0.7;
facecolor2 = [0.8500 0.3250 0.0980];
linecolor2 = facecolor2*0.7;

fig = figure('Position', [10 10 1600 1000]);
hold on
plot(translate(pgon1, [-1,1]), 'FaceColor', facecolor1, 'FaceAlpha', 0.5, 'LineStyle', '--', 'EdgeColor' , linecolor1);
plot(translate(pgon2, [0,0]), 'FaceColor', facecolor2, 'FaceAlpha', 0.5, 'LineStyle', '--', 'EdgeColor' , linecolor2);
hold off
plot_balls(t_opt, p_res, {1:2, 3:4}, [pgon1,pgon2], {facecolor1,facecolor2}, {linecolor1,linecolor2}, fig, 'tf_discs')

%%
if 1
    figure('Renderer', 'painters', 'Position', [100 100 1400 600])
    subplot(311)
    plot(t_opt,v1,'LineWidth',1.5);
    hold on
    plot(t_opt,v2,'LineWidth',1.5);
    legend({'$v_1(t)$','$v_2(t)$'},'interpreter','latex');
    xlabel('$t$','interpreter','latex');
    ylabel('$v(t)$','interpreter','latex');
    xlim([0,2])
    grid on
    % axis equal
    subplot(312)
    plot(t_opt,v3,'LineWidth',1.5);
    hold on
    plot(t_opt,v4,'LineWidth',1.5);
    grid on
    legend({'$v_3(t)$','$v_4(t)$'},'interpreter','latex');
    xlabel('$t$','interpreter','latex');
    ylabel('$v(t)$','interpreter','latex');
    xlim([0,2])
    subplot(313)
    stairs(t_opt(1:N_FE:end),[ocp_solver.get('u'),nan*ones(2,1)]','LineWidth',1.5);
    legend({'$u_1(t)$','$u_2(t)$'},'interpreter','latex');
    grid on
    xlabel('$t$','interpreter','latex');
    ylabel('$u$','interpreter','latex');
    xlim([0,2])
end