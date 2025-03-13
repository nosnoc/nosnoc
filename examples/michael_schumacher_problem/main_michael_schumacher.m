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

% This file is part of NOSNOC.%
clear; clc; close all
import casadi.*
%%
% An optimal control problem from:
% Optimal control of systems with discontinuous differential equations
% April 2012, Numerische Mathematik 114(4):653-695
% DOI: 10.1007/s00211-009-0262-2
% David E. Stewart Mihai Anitescu

%% Problem parameters
path_constraint = 'track';
track_width = 0.5;
omega = 1;
chicane_tightness = 1;
chicane_width = 2;

%% NOSNOC settings
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();

% Optimal control problem options
problem_options.n_s = 2;
problem_options.use_fesd = 1;
problem_options.cross_comp_mode = "FE_FE";
problem_options.T_final_max = 5*pi;
problem_options.T_final_min = 2;
problem_options.dcs_mode = 'Stewart';
problem_options.g_path_at_fe = 1;
problem_options.g_path_at_stg = 1;
problem_options.time_optimal_problem = 1;
problem_options.relax_terminal_constraint = ConstraintRelaxationMode.ELL_1;

% Solver options
solver_options.mpecopt.settings_casadi_nlp.ipopt.print_level = 5;
solver_options.sigma_0 = 1e0;
solver_options.opts_casadi_nlp.ipopt.max_iter = 1000;
solver_options.print_level = 5;
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';

%% Time discretization settings
problem_options.T = 5;
problem_options.N_stages = 30;
problem_options.N_finite_elements = 2;

%% Model equations
% Declare model variables
qx = SX.sym('qx'); qy = SX.sym('qy'); % position
vx = SX.sym('vx'); vy = SX.sym('vy'); % velocity
alpha = SX.sym('alpha'); % orientation
tangent = [cos(alpha);sin(alpha)];
normal = [-sin(alpha);cos(alpha)];

q = [qx;qy];
v = [vx;vy];
x = [q;v;alpha];
x0 = zeros(5,1);
%  control
a = SX.sym('a'); % acceleration
s = SX.sym('s'); % steering angle
u = [a;s];
u_max  = 2;
% Modes of the PSS
Friction_max = 4;
f_1 = [v;a*tangent+Friction_max*normal;s*(tangent'*v)];
f_2 = [v;a*tangent-Friction_max*normal;s*(tangent'*v)];


%% Populate nosnoc model object
model = nosnoc.model.Pss();
model.x0 = x0;
model.x = x;
model.lbx = [-10;-10;-inf;-inf;-inf];
model.ubx = [10;10;inf;inf;inf];
model.u = u;
model.lbu = -u_max*ones(2,1);
model.ubu = u_max*ones(2,1);
model.u0 = [0;0]; % inital guess for control variables
model.c = normal'*v; % switching function
model.S = [1;-1]; % f_1 is for c>0, f_2 is for c<0
model.F = [f_1 f_2];

%%  general nonlinear constraints for the track
switch path_constraint
    case 'none'
        q_target = [3*pi;0];
        xx = linspace(0,q_target(1),problem_options.N_stages);
        yy = 0*xx;
    case 'linear'
        model.g_path = qy-(qx);
        q_target = [3*pi;3*pi];
        xx = linspace(0,q_target(1),problem_options.N_stages);
        yy = xx;

    case 'nonlinear'
        model.g_path = qy-sin(omega*qx);
        q_target = [3*pi;sin(omega*3*pi)];
        xx = linspace(0,q_target(1),problem_options.N_stages);
        yy = sin(omega*xx);
    case 'track'
        arg1 = qx-pi;
        arg2 = qx-2*pi;
        sig = 1e-1;
        step1 = 0.5*(1+tanh(arg1/sig));
        step2 = 0.5*(1+tanh(arg2/sig));
        model.g_path = qy-(sin(qx)*(1-step1)+(pi-qx)*step1*(1-step2)+(-pi-sin(qx))*step2);
        q_target = [3*pi;-pi];
        xx = linspace(0,q_target(1),problem_options.N_stages);
        yy = sin(xx).*(1-0.5.*(1+tanh((xx-pi)/sig))) +...
            (pi-xx).*0.5.*(1+tanh((xx-pi)/sig)).*(1-0.5.*(1+tanh((xx-2*pi)/sig)))+...
            (-pi-sin(xx)).*0.5.*(1+tanh((xx-2.*pi)/sig));
    case 'chicane'

        q_target = [10;2*chicane_width];
        model.g_path = qy-((chicane_width)+chicane_width*tanh(chicane_tightness*(qx-q_target(1)/2)));
        xx = linspace(0,q_target(1),problem_options.N_stages);
        yy = (chicane_width)+chicane_width*tanh(chicane_tightness*(xx-q_target(1)/2));
end
if ~isequal(path_constraint,'none')
    model.lbg_path = -track_width;
    model.ubg_path = +track_width;
end


%% objective
if ~problem_options.time_optimal_problem
    model.f_q = u'*u;
    problem_options.T = 5;
else
    model.f_q = 0.01*u'*u;
end

% Terminal Constraint
model.g_terminal = q-q_target; % temrinal constraint expression

% intitial guess for fun
alpha = [atan2(diff(yy),diff(xx)),0];
x_guess = vertcat(xx,yy, zeros(2,problem_options.N_stages), alpha); % Not used ATM

%% Compute active set guess for mpecopt solver
active_set_guess = nosnoc.activeset.Pss({[1,2]},'times', [problem_options.T]); % [1,2] is a slinding mode (i.e. c = 0), and it is assumed to hold over [0,T]
solver_options.mpecopt.initialization_strategy = 'TakeProvidedActiveSet';

%% Solve
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.set_initial_active_set(active_set_guess); 
% [IG,IH,I00] = ocp_solver.discrete_time_problem.process_active_set(active_set_guess);
ocp_solver.solve('mpecopt');
fprintf('Objective values is: %2.4f \n', ocp_solver.get_objective());
%% Extract results
u_opt = ocp_solver.get("u");
u1_opt = u_opt(1,:);
u2_opt = u_opt(2,:);

x_opt = ocp_solver.get("x");
x1_opt = x_opt(1,:);
x2_opt = x_opt(2,:);
x3_opt = x_opt(3,:);
x4_opt = x_opt(4,:);
x5_opt = x_opt(5,:);
t_grid = ocp_solver.get_time_grid();
u_grid = ocp_solver.get_control_grid();

%% Plot all results
% controls
figure
stairs(u_grid,[u1_opt,nan])
hold on
stairs(u_grid,[u2_opt,nan])
xlabel('$t$','interpreter','latex');
ylabel('$u(t)$','interpreter','latex');
ylim([-2.2 2.2])
grid on

% plots
figure
subplot(211)
plot(x1_opt,x2_opt);
hold on
plot(x1_opt,x2_opt,'r.');
xlabel('$q_x$','interpreter','latex');
ylabel('$q_y$','interpreter','latex');
hold on
grid on
plot(q_target(1),q_target(2),'rx','MarkerSize',6)
xx = linspace(0,q_target(1),1e2);

switch path_constraint
    case 'none'
    case 'linear'
        plot(xx ,(xx)-track_width,'k');
        plot(xx ,(xx)+track_width,'k');
    case 'nonlinear'
        plot(xx ,sin(omega*xx)-track_width,'k');
        plot(xx ,sin(omega*xx)+track_width,'k');
    case 'track'
        arg1 = xx-pi;
        arg2 = xx-2*pi;
        sig = 1e-1;
        step1 = 0.5*(1+tanh(arg1/sig));
        step2 = 0.5*(1+tanh(arg2/sig));
        yy = sin(xx).*(1-step1)+(pi-xx).*step1.*(1-step2)+(-pi-sin(xx)).*step2;
        plot(xx ,yy-track_width,'k');
        hold on
        plot(xx ,yy+track_width,'k');
    case 'chicane'
        yy = (chicane_width)+chicane_width*tanh(chicane_tightness*(xx-q_target(1)/2));
        plot(xx ,yy-track_width,'k');
        hold on
        plot(xx ,yy+track_width,'k');
end
grid on
axis equal

subplot(212)
stairs(x1_opt(1:problem_options.N_finite_elements(1):end),[u1_opt,nan]);
hold on
stairs(x1_opt(1:problem_options.N_finite_elements(1):end),[u2_opt,nan]);
xlabel('$u(q_x)$','interpreter','latex');
ylabel('$q_x$','interpreter','latex');
legend({'$a(x)$','$s(x)$'},'interpreter','latex');
ylim([-2.2 2.2])
grid on

% Compute normal and tangetial velocity;
tangent = [cos(x5_opt);sin(x5_opt)];
normal = [-sin(x5_opt);cos(x5_opt)];
v_opt = [x3_opt;x4_opt];
v_normal = dot(normal,v_opt);
v_tangent = dot(tangent,v_opt);

figure
plot(t_grid,v_normal);
hold on;
plot(t_grid,v_tangent);
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
legend({'$n^\top v$','$t^\top v$'},'interpreter','latex');

%% Plot intermediate iterations
res_opt = ocp_solver.discrete_time_problem.w.res;
X_iters = ocp_solver.stats.iter.X_outer;
% X_iters = ocp_solver.stats.iter.X_all;
n_iters = size(X_iters,2);
figure
for ii = 1:n_iters
    res_ii = X_iters(:,ii);
    ocp_solver.discrete_time_problem.w.res = res_ii;
    ocp_solver.get('x')
    x1_opt = x_opt(1,:);
    x2_opt = x_opt(2,:);
    x3_opt = x_opt(3,:);
    x4_opt = x_opt(4,:);
    x5_opt = x_opt(5,:);
    tangent = [cos(x5_opt);sin(x5_opt)];
    normal = [-sin(x5_opt);cos(x5_opt)];
    v_opt = [x3_opt;x4_opt];
    v_normal = dot(normal,v_opt);
    v_tangent = dot(tangent,v_opt);
    subplot(121)
    plot(x1_opt,x2_opt);
    hold on;
    % plot(x1_opt,x2_opt,'r.');
    xlabel('$q_x$','interpreter','latex');
    ylabel('$q_y$','interpreter','latex');
    hold on
    grid on
    plot(q_target(1),q_target(2),'rx','MarkerSize',6)
    xx = linspace(0,q_target(1),1e2);

    switch path_constraint
        case 'none'
        case 'linear'
            plot(xx ,(xx)-track_width,'k');
            plot(xx ,(xx)+track_width,'k');
        case 'nonlinear'
            plot(xx ,sin(omega*xx)-track_width,'k');
            plot(xx ,sin(omega*xx)+track_width,'k');
        case 'track'
            arg1 = xx-pi;
            arg2 = xx-2*pi;
            sig = 1e-1;
            step1 = 0.5*(1+tanh(arg1/sig));
            step2 = 0.5*(1+tanh(arg2/sig));
            yy = sin(xx).*(1-step1)+(pi-xx).*step1.*(1-step2)+(-pi-sin(xx)).*step2;
            plot(xx ,yy-track_width,'k');
            hold on
            plot(xx ,yy+track_width,'k');
        case 'chicane'
            yy = (chicane_width)+chicane_width*tanh(chicane_tightness*(xx-q_target(1)/2));
            plot(xx ,yy-track_width,'k');
            hold on
            plot(xx ,yy+track_width,'k');
    end
    grid on
    axis equal
    subplot(122)
    plot(t_grid,v_normal);
    hold on
    plot(t_grid,v_tangent);
    xlabel('$t$','interpreter','latex');
    ylabel('$v$','interpreter','latex');
    legend({'$n^\top v$','$t^\top v$'},'interpreter','latex');
    grid on;
    pause(1)
end
ocp_solver.discrete_time_problem.w.res = res_opt; % restore solution
%% Plot via function
% plot_results_ms(model, problem_options, ocp_solver,...
% q_target, path_constraint, track_width, omega, chicane_tightness, chicane_width)
