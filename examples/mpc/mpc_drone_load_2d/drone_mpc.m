% Closed-loop MPC for the 2D drone-load example.
%
clear; close all;
import casadi.*
import nosnoc.*

%%  Flags
filename              = 'drone_load_mpc';
run_animation         = true;
use_ccopt             = 0;    % Recommended for best perfomance but needs CCOpt installation
use_rtmpc             = 1; 
save_video            = 0;
DT                    = 0.05;
T_pred                = 0.5;
T_sim                 = 4.0;
N_sim                 = 4;
model_plant_mismatch  = 0;

N_steps  = round(T_sim / DT);
N_stages = round(T_pred / DT);

%% Load drone model and user reference settings
[model, robot_data] = drone_load_model();

ref_cfg = struct();
% final horizontal target position of drone/load
ref_cfg.x_target = 6.0;
% carrying altitude of the drone
ref_cfg.z_drone_high = 3.5;
% nominal horizontal drone velocity during transfer
ref_cfg.vx_drone_ref = 1.0;
% final drone altitude after releasing the load
ref_cfg.z_drone_final = 0.55;

% fractions of T_sim used for climb / horizontal transport / descend
% must satisfy: climb + move < 1
ref_cfg.climb_fraction   = 0.1;
ref_cfg.move_fraction    = 0.75;
% descend phase is the remaining time

%% Reference schedule
u_ref_default = robot_data.hover_u;

%% Options
mpc_options = nosnoc.mpc.Options();
problem_options = nosnoc.Options();
homotopy_options = nosnoc.reg_homotopy.Options();

problem_options.N_stages = N_stages;
problem_options.N_finite_elements = ones(N_stages,1);
problem_options.h = DT;
problem_options.T = T_pred;
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 1;
problem_options.use_fesd = 0;
problem_options.friction_model = 'Polyhedral';
problem_options.ub_gamma_d = 1e3;
problem_options.ub_delta_d = 1e3;
problem_options.time_freezing = false;

homotopy_options.N_homotopy = 10;
homotopy_options.homotopy_update_slope = 0.1;
homotopy_options.complementarity_tol = 1e-5;
homotopy_options.sigma_0 = 1;
homotopy_options.relaxation_strategy = "SCHOLTES_INEQ";
homotopy_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
homotopy_options.opts_casadi_nlp.ipopt.max_iter = 1500;
homotopy_options.print_level = 3;

mpc_options.fast_sigma_0 = 1e-3;
mpc_options.solve_advanced_problem = 0;
mpc_options.do_shift_initialization = true;
mpc_options.warmstart_full_mpc = true;
mpc_options.warmstart_qpec = false;
mpc_options.use_probing_qp = false;
mpc_options.use_feedback_qp = false;
mpc_options.objective_ratio = 0.98;

%  CCOpt Solver Options
ccopt_options = nosnoc.ccopt.Options(); % ccopt options
ccopt_options.opts_madnlp.linear_solver = 'Ma27Solver';
% ccopt_options.opts_madnlp.tol = tol;
ccopt_options.opts_madnlp.barrier.TYPE = 'MonotoneUpdate';
ccopt_options.opts_madnlp.disable_garbage_collector = true;
% ccopt_options.opts_madnlp.barrier.mu_min = 1e-9;
% ccopt_options.opts_ccopt.relaxation_update.TYPE = 'RolloffRelaxationUpdate';
% ccopt_options.opts_ccopt.relaxation_update.rolloff_slope = 2.5;
% ccopt_options.opts_ccopt.relaxation_update.rolloff_point = 1e-12;
% ccopt_options.opts_ccopt.sigma_min = 1e-8;
%ccopt_options.opts_ccopt.relaxation_update.TYPE = 'RelaxLBUpdate';
%ccopt_options.opts_ccopt.relaxation_update.relax_threshold = 1e-5;
ccopt_options.opts_ccopt.q_regularization = 'critical_rho';
ccopt_options.opts_ccopt.critical_rho_factor = 0.9999;
%ccopt_options.opts_madnlp.barrier.mu_init = 1.0;
%ccopt_options.opts_madnlp.barrier.TYPE = 'QualityFunctionUpdate';

%% Model and parametric costs

if ~isfield(robot_data, 'x0')
    robot_data.x0 = full(model.x0);
end

n_x = length(model.x);
n_q = robot_data.n_q;
n_u = robot_data.n_u;
x = model.x;
u = model.u;

Qq = diag([20, 20, 2, 35, 35]);
Qv = diag([0.5, 0.5, 0.1, 0.5, 0.5]);
Q = blkdiag(Qq, Qv);
R = 5e-3 * eye(n_u);

Qq_T = 5*diag([80, 80, 10, 200, 100]);
Qv_T = diag([1, 1, 0.2, 1, 1]);
Q_T = blkdiag(Qq_T, Qv_T);

x_ref_sym = SX.sym('x_ref', n_x);
x_ref_T_sym = SX.sym('x_ref_T', n_x);
u_ref_sym = SX.sym('u_ref', n_u);

f_q = (x - x_ref_sym).' * Q * (x - x_ref_sym) + (u - u_ref_sym).' * R * (u - u_ref_sym);
f_q_T = (x - x_ref_T_sym).' * Q_T * (x - x_ref_T_sym);

model.f_q = f_q;
model.f_q_T = f_q_T;
model.p_global = [x_ref_T_sym; u_ref_sym]; % terminal state reference; constat control hover reference
model.p_time_var = x_ref_sym; % drone and load pos. and vel. are time-varying
model.p_global_val = [model.x0; u_ref_default];
model.p_time_var_val = repmat(model.x0, 1, N_stages);

%%  Set up MPC Solver (real-time or fully converged, ccopt or not)

if use_rtmpc
    if use_ccopt
        mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, ccopt_options, ccopt_options);
    else
        mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, homotopy_options, homotopy_options);
    end
    % mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, homotopy_options, mpecopt_options);
else
    if use_ccopt
        mpc = nosnoc.mpc.FullMpcc(model, mpc_options, problem_options, ccopt_options);
    else
        mpc = nosnoc.mpc.FullMpcc(model, mpc_options, problem_options, homotopy_options);
    end
end

%% Plant integrator
if model_plant_mismatch
    [sim_model, ~] = drone_load_model();
    sim_model.lbx = -inf(n_x,1);
    sim_model.ubx = inf(n_x,1);
    sim_model.lbu = model.lbu;
    sim_model.ubu = model.ubu;

    sim_problem_options = nosnoc.Options();
    integrator_options = nosnoc.integrator.Options();
    integrator_options.integrator_plugin = "FESD";

    sim_solver_options = integrator_options.fesd_solver_opts;
    sim_problem_options.rk_scheme = RKSchemes.RADAU_IIA;
    sim_problem_options.n_s = 1;
    sim_problem_options.N_finite_elements = 1;
    sim_problem_options.N_sim = N_sim;
    sim_problem_options.T_sim = DT;
    sim_problem_options.use_fesd = false;
    sim_problem_options.lift_velocity_state = 1;
    sim_problem_options.friction_model = 'Polyhedral';
    sim_problem_options.ub_gamma_d = 1e2;
    sim_problem_options.ub_delta_d = 1e2;

    sim_solver_options.print_level = 1;
    sim_solver_options.complementarity_tol = 1e-7;
    sim_solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';

    integrator = nosnoc.Integrator(sim_model, sim_problem_options, integrator_options);
end

%% Reference schedule
x_refs =  build_reference_table(robot_data, N_steps + N_stages + 1, DT, T_sim, ref_cfg);
u_ref_default = robot_data.hover_u;

for stage = 1:N_stages
    mpc.set_param('p_time_var', {stage}, x_refs(:, stage));
end
mpc.set_param('p_global', [], [x_refs(:, N_stages+1); u_ref_default]); % terminal reference and constant control ref

%% Closed loop
x0 = model.x0;
x_cl = x0;
t_cl = 0;
u_cl = [];
preparation_times = [];
feedback_times = [];

x_high_res = x0;
t_high_res = 0;
lambda_n_high_res = [];

for k = 1:N_steps
    % Online update of drone and load reference
    for stage = 1:N_stages
        mpc.set_param('p_time_var', {stage}, x_refs(:, k + stage));
    end
    mpc.set_param('p_global', [], [x_refs(:, k + N_stages + 1); u_ref_default]); % terminal reference

    stats_prep = mpc.do_preparation();
    preparation_times = [preparation_times, stats_prep.preparation_time]; %#ok<AGROW>

    [u_k, stats_fb] = mpc.get_feedback(x0);
    feedback_times = [feedback_times, stats_fb.feedback_time]; %#ok<AGROW>

    if model_plant_mismatch
        integrator.set_x0(x0);
        [t_grid_sim, x_sim] = integrator.simulate("u", repmat(u_k, [1, N_sim]), "x0", x0);
        x0 = x_sim(:, end);

        x_high_res = [x_high_res, x_sim(:,2:end)]; %#ok<AGROW>
        t_high_res = [t_high_res, t_high_res(end) + t_grid_sim(2:end)]; %#ok<AGROW>

        lambda_step = integrator.get('lambda_normal');
        lambda_n_high_res = [lambda_n_high_res, lambda_step]; %#ok<AGROW>
    else
        x0 = mpc.get_predicted_state();
    end

    x_cl = [x_cl, x0]; %#ok<AGROW>
    t_cl = [t_cl, t_cl(end) + DT]; %#ok<AGROW>
    u_cl = [u_cl, u_k]; %#ok<AGROW>

    fprintf('MPC step %d / %d\n', k, N_steps);
end

if model_plant_mismatch
    t_anim = t_high_res;
    x_anim = x_high_res;
else
    t_anim = t_cl;
    x_anim = x_cl;
end

q_anim = x_anim(1:n_q, :);
plot_limits = compute_plot_limits(q_anim, x_refs(1:n_q,:), robot_data, ref_cfg);

%% Plots
figure('Color','w');

subplot(3,1,1)
plot(t_anim, q_anim(1,:), 'LineWidth', 1.5); hold on
plot(t_anim, q_anim(2,:), 'LineWidth', 1.5);
plot(t_anim, q_anim(4,:), '--', 'LineWidth', 1.5);
plot(t_anim, q_anim(5,:), '--', 'LineWidth', 1.5);
legend({'x_d','z_d','x_l','z_l'}, 'Location', 'best');
ylabel('positions');
grid on

subplot(3,1,2)
plot(t_anim, q_anim(3,:), 'LineWidth', 1.5);
ylabel('\phi');
grid on

subplot(3,1,3)
stairs(t_cl(1:end-1), u_cl.', 'LineWidth', 1.2);
xlabel('t [s]');
ylabel('u');
grid on

if ~isempty(lambda_n_high_res)
    figure('Color','w');
    plot(lambda_n_high_res.', 'LineWidth', 1.2);
    legend({'\lambda_{tether}','\lambda_{ground}'}, 'Location', 'best');
    xlabel('index');
    ylabel('normal force');
    grid on
end

if run_animation || save_video
    animate_drone_load(t_anim, q_anim, robot_data, filename, save_video, plot_limits, ref_cfg);
end


%% Local functions
function x_refs = build_reference_table(robot_data, N_points, DT, T_sim, ref_cfg)
x0 = robot_data.x0;
n_q = robot_data.n_q;

load_z_ground = robot_data.ground_z + robot_data.load_h/2;

x_target = ref_cfg.x_target;
z_drone_high = max(ref_cfg.z_drone_high, load_z_ground + robot_data.L + 0.05);
z_load_high = z_drone_high - robot_data.L;
z_drone_final = ref_cfg.z_drone_final;

% phase times are defined on [0, T_sim], not on prediction horizon
t0 = 0;
t1 = ref_cfg.climb_fraction * T_sim;
t2 = (ref_cfg.climb_fraction + ref_cfg.move_fraction) * T_sim;
t3 = T_sim;

if t2 >= t3
    error('Need climb_fraction + move_fraction < 1.');
end

% optional consistency check against desired horizontal velocity
dx = abs(x_target - x0(1));
t_move_nom = dx / max(abs(ref_cfg.vx_drone_ref), 1e-6);
t_move_ref = t2 - t1;
fprintf('Reference move time = %.3f s, implied vx = %.3f m/s, requested vx ~= %.3f m/s\n', ...
    t_move_ref, dx / max(t_move_ref,1e-6), ref_cfg.vx_drone_ref);

% mission nodes on [0, T_sim]
t_nodes = [t0, t1, t2, t3];
q_nodes = [ ...
    x0(1), x0(1),         x_target,      x_target;       % drone x
    x0(2), z_drone_high,  z_drone_high,  z_drone_final; % drone z
    0.0,   0.0,           0.0,           0.0;           % drone pitch
    x0(4), x0(4),         x_target,      x_target;      % load x
    x0(5), z_load_high,   z_load_high,   load_z_ground  % load z
    ];

t_grid = (0:(N_points-1)) * DT;

% build ref only over [0, T_sim], then hold constant afterwards
q_refs = zeros(n_q, N_points);
for k = 1:N_points
    tk = min(t_grid(k), T_sim);
    q_refs(:,k) = interp1(t_nodes.', q_nodes.', tk, 'linear').';
end

v_refs = zeros(n_q, N_points);

% optional nonzero x velocity during transfer phase
vx_transfer = (x_target - x0(1)) / max(t2 - t1, 1e-6);
for k = 1:N_points
    tk = min(t_grid(k), T_sim);
    if tk >= t1 && tk <= t2
        v_refs(1,k) = vx_transfer; % drone vx
        v_refs(4,k) = vx_transfer; % load vx
    end
end

% after T_sim hold the landed final state with zero velocity
idx_hold = t_grid > T_sim;
if any(idx_hold)
    q_refs(:, idx_hold) = repmat(q_nodes(:,end), 1, nnz(idx_hold));
    v_refs(:, idx_hold) = 0;
end

x_refs = [q_refs; v_refs];
end

function plot_limits = compute_plot_limits(q_traj, q_refs, robot_data, ref_cfg)
x_all = [q_traj(1,:), q_traj(4,:), q_refs(1,:), q_refs(4,:), ref_cfg.x_target];
z_all = [q_traj(2,:), q_traj(5,:), q_refs(2,:), q_refs(5,:), robot_data.ground_z];

xmin_data = min(x_all) - 0.5 * max(robot_data.drone_w, robot_data.load_w);
xmax_data = max(x_all) + 0.5 * max(robot_data.drone_w, robot_data.load_w);

zmin_data = min([robot_data.ground_z, z_all - 0.5 * robot_data.load_h]);
zmax_data = max(z_all) + 0.5 * robot_data.drone_h;

span_x = max(xmax_data - xmin_data, 1.0);
span_z = max(zmax_data - zmin_data, 1.0);

pad_x = max(0.2, 0.12 * span_x);
pad_z = max(0.15, 0.15 * span_z);

xmin = xmin_data - pad_x;
xmax = xmax_data + pad_x;
zmin = min(robot_data.ground_z - 0.05, zmin_data - 0.05);
zmax = zmax_data + pad_z;

plot_limits = [xmin, xmax, zmin, zmax];
end

function animate_drone_load(t, q_traj, robot_data, filename, save_video, plot_limits, ref_cfg)
N = size(q_traj, 2);

fig = figure('Color', 'w');

if save_video
    video = VideoWriter([filename '.mp4'], 'MPEG-4');
    open(video);
end

for k = 1:N
    clf(fig);
    hold on
    axis equal
    xlim(plot_limits(1:2));
    ylim(plot_limits(3:4));
    xlabel('x');
    ylabel('z');

    plot([plot_limits(1)-1, plot_limits(2)+1], ...
        [robot_data.ground_z, robot_data.ground_z], ...
        'k-', 'LineWidth', 1.5);

    plot([ref_cfg.x_target, ref_cfg.x_target], ...
        [robot_data.ground_z, min(plot_limits(4), robot_data.ground_z + 0.35)], ...
        'k--', 'LineWidth', 1.0);

    p_d = q_traj(1:2, k);
    phi = q_traj(3, k);
    p_l = q_traj(4:5, k);

    d = norm(p_d - p_l);
    if d >= robot_data.L - 1e-3
        plot([p_d(1), p_l(1)], [p_d(2), p_l(2)], 'k-', 'LineWidth', 2.0);
    else
        plot([p_d(1), p_l(1)], [p_d(2), p_l(2)], 'k--', 'LineWidth', 1.5);
    end

    draw_drone(p_d, phi, robot_data);
    draw_load(p_l, robot_data);

    title(sprintf('t = %.2f s', t(k)));
    drawnow

    if save_video
        writeVideo(video, getframe(fig));
    end
end

if save_video
    close(video);
end
end

function draw_drone(p_d, phi, robot_data)
R = [cos(phi), -sin(phi); sin(phi), cos(phi)];

body = 0.5 * [ -robot_data.drone_w,  robot_data.drone_w,  robot_data.drone_w, -robot_data.drone_w;
    -robot_data.drone_h, -robot_data.drone_h,  robot_data.drone_h,  robot_data.drone_h ];
body_w = R * body + p_d;

fill(body_w(1,:), body_w(2,:), [0.2, 0.45, 0.9], ...
    'EdgeColor', 'k', 'LineWidth', 1.2);

rotor_x = 0.5 * [-robot_data.drone_w, -robot_data.drone_w/3, robot_data.drone_w/3, robot_data.drone_w];
rotor_z = 0.55 * robot_data.drone_h * ones(1,4);
rotors = R * [rotor_x; rotor_z] + p_d;

plot(rotors(1,:), rotors(2,:), 'ko', ...
    'MarkerFaceColor', 'k', 'MarkerSize', 4);
end

function draw_load(p_l, robot_data)
rect = 0.5 * [ -robot_data.load_w,  robot_data.load_w,  robot_data.load_w, -robot_data.load_w;
    -robot_data.load_h, -robot_data.load_h,  robot_data.load_h,  robot_data.load_h ];
rect_w = rect + p_l;

fill(rect_w(1,:), rect_w(2,:), [0.85, 0.45, 0.2], ...
    'EdgeColor', 'k', 'LineWidth', 1.2);
end