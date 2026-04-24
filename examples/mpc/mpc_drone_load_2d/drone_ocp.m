%DRONE_OCP  Open-loop OCP for a 2D drone transporting a load with a tether.
%
%   The scenario is:
%   - initial: load on the ground, tether taut, drone low,
%   - middle: drone flies higher and transports the load,
%   - final: load is placed at a target position on the ground, drone descends
%     so the tether is slack (distance < L).
close all; clear; 
import casadi.*
import nosnoc.*

filename       = 'drone_load_ocp';
run_animation = true;
save_video    = false;


%% NOSNOC options
problem_options = nosnoc.Options();
solver_options  = nosnoc.reg_homotopy.Options();

solver_options.N_homotopy = 6;
solver_options.homotopy_update_slope = 0.1;
solver_options.complementarity_tol = 1e-6;
solver_options.sigma_0 = 1;
solver_options.relaxation_strategy = "SCHOLTES_INEQ";
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
solver_options.opts_casadi_nlp.ipopt.max_iter = 1e3;
solver_options.print_level = 3;

problem_options.rk_scheme = RKSchemes.RADAU_IIA;

problem_options.cross_comp_mode = 3;
problem_options.lift_velocity_state = 0;
problem_options.friction_model = 'Polyhedral';
problem_options.ub_gamma_d = 1e3;
problem_options.ub_delta_d = 1e3;
problem_options.time_freezing = 0;

if problem_options.time_freezing
    dt            = 0.2;
    N_stages      = 25;
    problem_options.use_fesd = 1;
    problem_options.n_s = 2;
    problem_options.N_finite_elements = 3;
    
else
    dt            = 0.05;
    N_stages      = 100;
    problem_options.use_fesd = 0; % If true, uses FESD-J; harder problem but much higher accuracy due to switch detection (consider also larger N_FE and n_s)
    problem_options.n_s = 1;
    problem_options.N_finite_elements = 1;
end

problem_options.T = N_stages * dt;
problem_options.N_stages = N_stages;
   

%% Model
[model, robot_data] = drone_load_model();
if problem_options.time_freezing
    model.dims.n_dim_contact = 1; % Time-freezing needs to know about the tangent space dimension
end
n_q = robot_data.n_q;
n_u = robot_data.n_u;
x = model.x;
u = model.u;

[q_ref, x_ref, x_ref_T] = build_reference(robot_data, problem_options.T, N_stages);

Qq = diag([20, 20, 2, 35, 40]);
Qv = diag([0.5, 0.5, 0.1, 0.5, 0.5]);
Q = blkdiag(Qq, Qv);
R = 5e-3 * eye(n_u);
Qq_T = diag([80, 80, 10, 150, 200]);
Qv_T = diag([1, 1, 0.2, 1, 1]);
Q_T = blkdiag(Qq_T, Qv_T);

u_ref = robot_data.hover_u;
model.lsq_x = {x, x_ref, Q};
model.lsq_u = {u, u_ref, R};
model.lsq_T = {x, x_ref_T, Q_T};

%% Solve
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();

x_res = ocp_solver.get('x');
u_res = ocp_solver.get('u');
if ~problem_options.time_freezing
    lambda_n = ocp_solver.get('lambda_normal');
end

q_res = x_res(1:n_q,:);
v_res = x_res(n_q+1:end,:);
t_x = linspace(0, problem_options.T, size(x_res,2));
t_u = linspace(0, problem_options.T - dt, size(u_res,2));

%% Plots
figure('Color','w');
subplot(3,1,1)
plot(t_x, q_res(1,:), 'LineWidth', 1.5); hold on
plot(t_x, q_res(2,:), 'LineWidth', 1.5);
plot(t_x, q_res(4,:), '--', 'LineWidth', 1.5);
plot(t_x, q_res(5,:), '--', 'LineWidth', 1.5);
legend({'x_d','z_d','x_l','z_l'}, 'Location', 'best');
ylabel('positions'); grid on

subplot(3,1,2)
plot(t_x, q_res(3,:), 'LineWidth', 1.5);
ylabel('\phi'); grid on

subplot(3,1,3)
stairs(t_u, u_res.', 'LineWidth', 1.2);
xlabel('t [s]'); ylabel('u'); grid on


if ~problem_options.time_freezing
    figure('Color','w');
    plot(lambda_n.', 'LineWidth', 1.5);
    legend({'\lambda_{tether}','\lambda_{ground}'}, 'Location', 'best');
    xlabel('index');
    ylabel('normal force');
    grid on
end

if run_animation || save_video
    animate_drone_load(t_x, q_res, robot_data, filename, save_video);
end

%% Helpers


function [q_ref, x_ref, x_ref_T] = build_reference(robot_data, T, N_stages)
x0 = robot_data.x0;
load_z_ground = robot_data.ground_z + robot_data.load_h/2;
x_target = 4.0;
z_drone_high = 2.55;
z_load_high = z_drone_high - robot_data.L;
z_drone_final = 0.45;

t_nodes = [0, 0.25*T, 0.75*T, T];
q_nodes = [x0(1), x0(1),     x_target, x_target;
    x0(2), z_drone_high, z_drone_high, z_drone_final;
    0.0,   0.0,       0.0,      0.0;
    x0(4), x0(4),     x_target, x_target;
    x0(5), z_load_high, z_load_high, load_z_ground];

t_grid = linspace(0, T, N_stages);
q_ref = interp1(t_nodes.', q_nodes.', t_grid, 'linear').';
v_ref = zeros(robot_data.n_q, N_stages);
x_ref = [q_ref; v_ref];
x_ref_T = [q_nodes(:,end); zeros(robot_data.n_q,1)];
end

function animate_drone_load(t, q_traj, robot_data, filename, save_video)
N = size(q_traj, 2);

fig = figure('Color', 'w');
axis equal
hold on
xlabel('x'); ylabel('z');
xlim([-0.6, 5.2]);
ylim([-0.05, 3]);

if save_video
    video = VideoWriter([filename '.mp4'], 'MPEG-4');
    if numel(t) > 1
        video.FrameRate = max(5, round(1 / max(t(2)-t(1), 1e-3)));
    else
        video.FrameRate = 20;
    end
    open(video);
end

for k = 1:N
    cla
    hold on
    axis equal
    xlim([-0.6, 5.2]);
    ylim([-0.05, 3]);

    plot([-1, 6], [robot_data.ground_z, robot_data.ground_z], 'k-', 'LineWidth', 1.5);
    plot([4.0, 4.0], [robot_data.ground_z, 0.35], 'k--', 'LineWidth', 1.0);

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
fill(body_w(1,:), body_w(2,:), [0.2, 0.45, 0.9], 'EdgeColor', 'k', 'LineWidth', 1.2);

rotor_x = 0.5 * [-robot_data.drone_w, -robot_data.drone_w/3, robot_data.drone_w/3, robot_data.drone_w];
rotor_z = 0.55 * robot_data.drone_h * ones(1,4);
rotors = R * [rotor_x; rotor_z] + p_d;
plot(rotors(1,:), rotors(2,:), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
end

function draw_load(p_l, robot_data)
rect = 0.5 * [ -robot_data.load_w,  robot_data.load_w,  robot_data.load_w, -robot_data.load_w;
    -robot_data.load_h, -robot_data.load_h,  robot_data.load_h,  robot_data.load_h ];
rect_w = rect + p_l;
fill(rect_w(1,:), rect_w(2,:), [0.85, 0.45, 0.2], 'EdgeColor', 'k', 'LineWidth', 1.2);
end

function value = get_opt(s, name, default)
if isstruct(s) && isfield(s, name)
    value = s.(name);
else
    value = default;
end
end
