%% Spinner OCP
clear; clc; close all;
import casadi.*

filename = 'spinner_traj';
run_animation = 1;

%% Default NOSNOC settings
problem_options = nosnoc.Options();
solver_options  = nosnoc.reg_homotopy.Options();

solver_options.N_homotopy = 6;
solver_options.homotopy_update_slope = 0.1;
solver_options.complementarity_tol = 1e-5;
solver_options.sigma_0 = 1;
solver_options.relaxation_strategy = "SCHOLTES_INEQ";
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
solver_options.opts_casadi_nlp.ipopt.max_iter = 800;
solver_options.print_level = 3;

problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 1;
problem_options.cross_comp_mode = 3;
problem_options.use_fesd = 0;
problem_options.lift_velocity_state = 1;
problem_options.friction_model = 'Polyhedral';

problem_options.ub_gamma_d = 1e3;
problem_options.ub_delta_d = 1e3;

%% Load spinner model

shape = 'ellipse';  % options: 'circle' (default) or 'ellipse'
[model, robot_data] = spinner_model(shape);
% model.mu = 1.2;

% model.J_tangent = [model.J_tangent, -model.J_tangent];
model.D_tangent = [model.J_tangent, -model.J_tangent];
n_q = robot_data.n_q; n_u = robot_data.n_u;
x0 = model.x0;

%% Initial and target states
q_init = [0.2; 1.5; 0.0];
v_init = [0;0;0];
x0 = [q_init; v_init];

q_nom_start = [0.3; 1.5; 0.0];
q_nom_end   = [0.3; 1.5; 6*pi];
v_nom_end   = [0;0;0];

%% Time horizon and discretization
dt = 0.05;
N_stg = 200;
T = N_stg * dt;
problem_options.T = T;
problem_options.N_stages = N_stg;
problem_options.N_finite_elements = 1;

%% Reference trajectory
q_ref = interp1([0 1],[q_nom_start q_nom_end]',linspace(0,1,N_stg),'linear')';
v_ref = repmat(v_nom_end,1,N_stg);
x_ref = [q_ref; v_ref];

% x_ref = [q_nom_end; v_nom_end];

%% Weights
Qq = diag([0.01 0.01 0.01]);
Qv = diag([0.1 0.1 0.1]);
Q  = blkdiag(Qq, Qv);
R  = diag([0.1 0.1])*1e-1;
Qfq = diag([10 10 10]);
Qfv = diag([0.1 0.1 0.1]);
Q_terminal = blkdiag(Qfq, Qfv);

u_ref = zeros(n_u,1);

%% Populate model cost terms
x = model.x; u = model.u;
model.lsq_x = {x, x_ref, Q};
model.lsq_u = {u, u_ref, R};
model.lsq_T = {x, [q_nom_end; v_nom_end], Q_terminal};

%% Call NOSNOC OCP solver
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();

%% Retrieve results
x = ocp_solver.get('x');
u = ocp_solver.get('u');
q_res = x(1:n_q,:);
v_res = x(n_q+1:end,:);
lambda_n = ocp_solver.get('lambda_normal');
lambda_t = ocp_solver.get('lambda_tangent');

%% Visualization
if run_animation
    h = figure('Color','w');
    hold on; axis equal
    xlabel('$x$','Interpreter','latex');
    ylabel('$z$','Interpreter','latex');
    xlim([-0.2 2.2]); ylim([-0.2 2.0]);

    % Extract geometry
    l1 = robot_data.l1;
    l2 = robot_data.l2;
    R  = robot_data.R;
    cxy = robot_data.cxy;
    N_frames = size(q_res,2);
    cx = robot_data.cxy(1);
    cz = robot_data.cxy(2);

    outputVideo = VideoWriter([filename '.mp4'],'MPEG-4');
    outputVideo.FrameRate = 1/dt/1; % playback slower
    open(outputVideo);

    for k = 1:N_frames
        cla
        q1 = q_res(1,k); q2 = q_res(2,k); q3 = q_res(3,k);

        % Joint positions
        p0 = [0;0];
        p1 = p0 + [l1*sin(q1); l1*cos(q1)];
        p2 = p1 + [l2*sin(q1+q2); l2*cos(q1+q2)];
        p_ball = cxy;

        % % Draw spinner (ball)
        % theta = linspace(0,2*pi,40);
        % plot(p_ball(1)+R*cos(theta), p_ball(2)+R*sin(theta), 'b','LineWidth',2);
        % % Ball rotation marker
        % m_ang = q3;
        % plot(p_ball(1)+R*cos(m_ang), p_ball(2)+R*sin(m_ang), 'go','MarkerFaceColor','g');

        switch lower(robot_data.shape)
            case 'circle'
                theta = linspace(0,2*pi,60);
                R = robot_data.a;     % same since a=b=R
                plot(cx + R*cos(theta), cz + R*sin(theta), 'b','LineWidth',2);
                m_ang = q3;
                plot(cx + R*cos(m_ang), cz + R*sin(m_ang), 'go','MarkerFaceColor','g');

            case 'ellipse'
                a  = robot_data.a;      b  = robot_data.b;
                q3 = q_res(3,k);  % or x_input(3) if you use x_input

                th = linspace(0, 2*pi, 200);                 % 1×N
                Xloc = a*cos(th);                            % 1×N
                Zloc = b*sin(th);                            % 1×N

                R2 = [cos(q3) -sin(q3);                      % 2×2 rotation
                    sin(q3)  cos(q3)];

                P = R2 * [Xloc; Zloc];                       % 2×N
                Xw = P(1,:) + cx;                            % 1×N
                Zw = P(2,:) + cz;                            % 1×N

                plot(Xw, Zw, 'b', 'LineWidth', 2); hold on

                % orientation marker: local +x end rotated to world
                Xmark = R2 * [a; 0];
                % plot(cx + Xmark(1), cz + Xmark(2), 'go', 'MarkerFaceColor', 'g');
        end


        % Draw arm
        plot([p0(1) p1(1)], [p0(2) p1(2)], 'r','LineWidth',6);
        plot([p1(1) p2(1)], [p1(2) p2(2)], 'r','LineWidth',6);
        plot(p0(1), p0(2), 'ko','MarkerFaceColor','k');
        plot(p1(1), p1(2), 'ko','MarkerFaceColor','k');
        plot(p2(1), p2(2), 'ko','MarkerFaceColor','k');
        plot(p_ball(1), p_ball(2), 'bo','MarkerFaceColor','b');

        title(sprintf('t = %.2f s', (k-1)*dt));
        drawnow
        frame = getframe(h);
        writeVideo(outputVideo, frame);
    end
    close(outputVideo);
end

%% figure
figure
subplot(121)
plot(q_res(3,:))
subplot(122)
plot(lambda_n)
hold on
plot(lambda_t')
