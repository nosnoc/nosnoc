clear all; clc; close all;

import casadi.*
import nosnoc.*

%% Problem size
N_blocks = 7;   % increase the number of blocks to have more states (to scale up the problem)

%% Time discretization / FESD options
N_FE    = 3;
use_fesd = 1;
gamma_h = 1.0;

dt       = 0.1;      % time step
T        = 10.0;     % total horizon
T        = 5.0;     % total horizon
N_stages = T/dt;

use_f0 = true;

%% Model parameters
mu_fric  = 0.4;      % Coulomb friction magnitude (same for all blocks)
k_spring = 5.0;      % spring stiffness
d_damp   = 0.1;      % viscous damping on velocities

%% Equilibrium positions
% q_eq(i) = (i-1)*0.2
q_eq = (0:(N_blocks-1))' * 0.2;

%% States: positions and velocities of all blocks
q = SX.sym('q', N_blocks, 1);     % absolute positions of blocks
v = SX.sym('v', N_blocks, 1);     % velocities
x = [q; v];
n_x = length(x);

%% Control: force on first block only
u = SX.sym('u');    % scalar control

%% Build smooth part of dynamics (without Coulomb friction)
qdot = v;

% vdot from spring-damper chain with:
% - spring between wall at 0 and block 1
% - springs between neighboring blocks
F_spring = SX.zeros(N_blocks,1);

% Wall–block 1 spring
s_wall = (q(1) - 0) - (q_eq(1) - 0);           % displacement from rest length
F_spring(1) = F_spring(1) - k_spring * s_wall;

% Springs between blocks i-1 and i
for i = 2:N_blocks
    s_ij = (q(i) - q(i-1)) - (q_eq(i) - q_eq(i-1));  % displacement from rest
    F_spring(i)   = F_spring(i)   - k_spring * s_ij; % force on block i
    F_spring(i-1) = F_spring(i-1) + k_spring * s_ij; % equal and opposite on i-1
end

% Add damping and control
vdot = F_spring - d_damp * v;
vdot(1) = vdot(1) + u;     % control acts only on block 1

% Base dynamics (without Coulomb friction)
f_base = [qdot; vdot];

%% Build nonsmooth Coulomb friction using PSS structure
% Switching functions: each block's velocity
c_cells = cell(1, N_blocks);
S_cells = cell(1, N_blocks);
F_cells = cell(1, N_blocks);

for i = 1:N_blocks
    % switching function c_i = v_i
    c_cells{i} = v(i);
    S_cells{i} = [1; -1];

    % index in the state derivative corresponding to vdot_i
    idx_vdot_i = N_blocks + i;

    % friction contribution: +/- mu_fric in vdot_i
    d_minus = SX.zeros(n_x,1);
    d_plus  = SX.zeros(n_x,1);

    d_minus(idx_vdot_i) = -mu_fric;
    d_plus(idx_vdot_i)  =  mu_fric;

    if i == 1
        % first "layer" carries the base dynamics
        if ~use_f0 
            f_i1 = f_base + d_minus;
            f_i2 = f_base + d_plus;
        else
            f_i1 = d_minus;
            f_i2 = d_plus;
        end
    else
        % other layers only contribute the friction offsets
        f_i1 = d_minus;
        f_i2 = d_plus;
    end

    F_cells{i} = [f_i1 f_i2];
end



%% Bounds and initial condition
% Initial positions: slightly perturbed from equilibrium
q0 = q_eq + 0.2;
v0 = zeros(N_blocks,1);
x0 = [q0; v0];

% State bounds (wide enough around eq)
q_min = min(q_eq) - 2.0;
q_max = max(q_eq) + 2.0;
v_max = 10.0;
lbx = [q_min*ones(N_blocks,1); -v_max*ones(N_blocks,1)];
ubx = [q_max*ones(N_blocks,1);  v_max*ones(N_blocks,1)];

% Control bounds
u_max = 20.0;
lbu = -u_max;
ubu =  u_max;

%% References and costs
% Reference = equilibrium positions and zero velocities
q_ref = q_eq;
v_ref = zeros(N_blocks,1);
x_ref = [q_ref; v_ref];
u_ref = 0;

% Weights: track positions to eq, then velocities, then control
Q_diag = [10*ones(N_blocks,1); 0.1*ones(N_blocks,1)];
Q_cost = diag(Q_diag);
R_cost = 0.01;
P_cost = 10 * Q_cost;   % terminal weight

%% NOSNOC model
model = nosnoc.model.Pss();
model.x   = x;
model.x0  = x0;
model.lbx = lbx;
model.ubx = ubx;

model.u   = u;
model.lbu = lbu;
model.ubu = ubu;

model.S = S_cells;
model.c = c_cells;
model.F = F_cells;

if use_f0 
    model.f_0 = f_base;
end

% Stage and terminal cost
model.f_q   = (x - x_ref).' * Q_cost * (x - x_ref) + (u - u_ref).' * R_cost * (u - u_ref);
model.f_q_T = (x - x_ref).' * P_cost * (x - x_ref);

%% NOSNOC options (OCP)
problem_options = nosnoc.Options();
solver_options  = nosnoc.reg_homotopy.Options();
ccopt_options  = nosnoc.ccopt.Options();

% RK / FESD settings
problem_options.rk_scheme          = RKSchemes.RADAU_IIA;
problem_options.n_s                = 2;
problem_options.cross_comp_mode    = "FE_FE";
problem_options.N_stages           = N_stages;
problem_options.N_finite_elements  = N_FE;
problem_options.T                  = T;
problem_options.use_fesd           = use_fesd;
problem_options.gamma_h            = gamma_h;

% Solver options
solver_options.complementarity_tol = 1e-7;
solver_options.N_homotopy          = 10;
solver_options.homotopy_steering_strategy ="DIRECT";
solver_options.lift_complementarities = 0;
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'mumps';

ccopt_options.opts_madnlp.tol = 1e-7;
ccopt_options.opts_madnlp.linear_solver = 'Ma27Solver';
ccopt_options.opts_ccopt.relaxation_update.TYPE = 'RolloffRelaxationUpdate';
ccopt_options.opts_ccopt.relaxation_update.rolloff_slope = 2.0;
ccopt_options.opts_ccopt.relaxation_update.rolloff_point = 1e-6;
ccopt_options.opts_ccopt.use_magic_step = false;

%% Create and solve OCP
%ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver = nosnoc.ocp.Solver(model, problem_options, ccopt_options);
ocp_solver.solve();
try
    wall_time = sum(ocp_solver.stats.wall_time)
    cpu_time = sum(ocp_solver.stats.cpu_time)
catch
    wall_time = ocp_solver.stats.t_wall_total
    cpu_time = ocp_solver.stats.t_proc_total
end
%% Extract solution
x_opt = ocp_solver.get('x');
u_opt = ocp_solver.get('u');
t_grid = ocp_solver.get_time_grid();
t_u    = ocp_solver.get_control_grid();

q_opt = x_opt(1:N_blocks,:);
v_opt = x_opt(N_blocks+1:end,:);

%% Plots: trajectories
figure

subplot(3,1,1)
plot(t_grid, q_opt', 'LineWidth', 1.5)
hold on;
plot(t_grid, repmat(q_eq',length(t_grid),1),'k--','LineWidth',1); % equilibrium lines
xlabel('$t$','Interpreter','latex')
ylabel('$q_i$','Interpreter','latex')
title('Positions')
grid on

subplot(3,1,2)
plot(t_grid, v_opt', 'LineWidth', 1.5)
hold on;
yline(0, 'k--', 'LineWidth', 1.0);
xlabel('$t$','Interpreter','latex')
ylabel('$v_i$','Interpreter','latex')
title('Velocities')
grid on

subplot(3,1,3)
stairs(t_u, [u_opt, u_opt(:,end)]', 'LineWidth', 1.5)
xlabel('$t$','Interpreter','latex')
ylabel('$u$','Interpreter','latex')
title('Control (first block)')
grid on

%% Simple 1D animation of blocks and springs
figure;
x_min = q_min - 0.5;
x_max = q_max + 0.5;
y_min = -0.5;
y_max =  0.5;

for k = 1:size(q_opt,2)
    clf; hold on; box on;

    % Wall at x = 0
    % plot([0 0],[y_min y_max],'k','LineWidth',2);
    % text(0,y_max,'Wall','HorizontalAlignment','left','VerticalAlignment','top');

    % Plot springs and blocks
    prev_x = 0;  % wall position
    for i = 1:N_blocks
        curr_x = q_opt(i,k);

        % Spring as a simple line between prev_x and curr_x
        plot([prev_x curr_x],[0 0],'LineWidth',1.5);

        % Block as a square marker
        plot(curr_x,0,'s','MarkerSize',12,'LineWidth',2);

        prev_x = curr_x;
    end

    % Equilibrium positions as small markers
    plot(q_eq, 0*q_eq + 0.1, 'ko', 'MarkerSize',5);

    xlim([x_min x_max]);
    ylim([y_min y_max]);
    title(sprintf('Time t = %.2f', t_grid(k)));
    xlabel('Position');
    drawnow;
end
