function [model,results] = run_quadroped_half_exp(scenario)

import casadi.*

% delete old gif
filename = scenario.filename;
run_animation = 1;
video_speed_up = 1;
objective = scenario.objective; % jump % wlk
%%
settings = NosnocOptions();
settings.irk_scheme = scenario.IrkScheme;
settings.dcs_mode = 'CLS';
settings.n_s = 2;  % number of stages in IRK methods
settings.N_homotopy = scenario.N_homotopy;
settings.opts_casadi_nlp.ipopt.tol = 1e-6;
settings.opts_casadi_nlp.ipopt.acceptable_tol = 1e-6;
settings.opts_casadi_nlp.ipopt.acceptable_iter = 3;

settings.opts_casadi_nlp.ipopt.max_iter = scenario.max_iter;
settings.equidistant_control_grid = 1;
settings.sigma_0 = scenario.sigma_0;
settings.cross_comp_mode = scenario.cross_comp_mode;
%% IF HLS solvers for Ipopt installed (check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions) use the settings below for better perfmonace:
settings.opts_casadi_nlp.ipopt.linear_solver = 'ma57';

%% discretization
model = half_unitri_ai_model();
unfold_struct(model,'caller');
if isequal(objective,'walk')
    T = scenario.T; % prediction horizon
    N_stg = scenario.N_stg; % control intervals
    N_FE = scenario.N_FE;  % integration steps per control intevral
    x_mid = x0;
    x_mid(2) = 0.4;
    x_end = x0;
    x_end(1) = scenario.q_x_final;
    x_mid(1) = (x_end(1)+x0(1))/2;
    Q = scenario.Q;
    Q_terminal= scenario.Q_terminal;
elseif isequal(objective,'jump')
    T = scenario.T; % prediction horizon
    N_stg = scenario.N_stg; % control intervals
    N_FE = scenario.N_FE;  % integration steps per control intevral
    x_mid = x0;
    x_end = x0;
    x_end(2) = scenario.q_z_final;
    x_mid(2) = (x_end(2)+x0(2))/2;
    Q = scenario.Q;
    Q_terminal= scenario.Q_terminal;
    if scenario.impose_terminal_constraint
        model.g_terminal = x(1:n_q)-x_end(1:n_q);
    end
else
end

u_ref = zeros(n_u,1);
R = scenario.R;
model.ubx(1) = scenario.q_x_final+1;
%% interpolate refernece
x_ref = interp1([0 0.5 1],[x0,x_mid,x_end]',linspace(0,1,N_stg),'spline')'; %spline

%% Populate model
model.mu = 0.8;
model.T = T;
model.dims.N_stages = N_stg;
model.dims.N_finite_elements  = N_FE;
% LSQ objective
model.lsq_x = {x,x_ref,Q};
model.lsq_u = {u,u_ref,R};
model.lsq_T = {x,x_end,Q_terminal};
model.u0 = 0*ones(n_u,1);
%% Call nosnoc solver
solver = NosnocSolver(model, settings);
[results,stats] = solver.solve();
%% read and plot results
unfold_struct(results,'caller');
%
q_res = x_opt(1:n_q,:);
v_res = x_opt(n_q+1:2*n_q,:);
N_frames = length(q_res);
h = T/(N_stg*N_FE);


if run_animation
    ground_X = [-0.1 2];
    ground_Y = [0  0];
    Quadruped_X = zeros(8, N_frames);
    Quadruped_Y = zeros(8, N_frames);
    for ii = 1:N_frames
        if ii == 1
            x_input = x0;
        else
            x_input = x_opt(:, ii);
        end
        [Quadruped_X(:, ii), Quadruped_Y(:, ii)] = getHalfQuadrupedConfiguration(model, x_input(1:11, 1));
    end

    %%
    ground_X = [min(q_res(1,:))-0.5  max(q_res(1,:))+0.5];
    ground_patch_X = [ground_X  ground_X(end:-1:1)];
    ground_patch_Y = [ground_Y-2  ground_Y];
    torso_color = [0 0 0];
    back_leg_color = 0.6*[1 1 1];
    front_leg_color = 0.3*[1 1 1];
    torso_line_width = 8;
    leg_up_line_width = 5;
    leg_down_line_width = 2;
    ground_brightness = 0.8;

    figure('Renderer', 'painters', 'Position', [10 10 900 600])
    pos = get(gcf, 'Position');
    % on 4k
    width = pos(3)*1.5;
    height = pos(4)*1.5;

%     width = pos(3)*1.0;
%     height = pos(4)*1.0;
    frames_gif = zeros(height, width, 1, N_frames, 'uint8');
    frames_video = { };
    for ii = 1:N_frames
        base = plot(ground_X, ground_Y, '.-k', 'MarkerSize', 1, 'LineWidth', 2);
        axis equal
        xlim(ground_X);
        ylim([-0.2 1.5])
        xlabel('$x$','Interpreter','latex');
        ylabel('$z$','Interpreter','latex');
        patch(ground_patch_X,ground_patch_Y,ground_brightness*ones(1,3));
        % torso
        hold on
        Quadruped_torso = plot(Quadruped_X(1:2, ii), Quadruped_Y(1:2, ii), '.-','Color', torso_color, 'MarkerSize', 15, 'LineWidth', torso_line_width);% torso
        % fron legs
        hold on
        Quadruped_leg1_up = plot(Quadruped_X(3:4, ii), Quadruped_Y(3:4, ii), '.-', 'Color', front_leg_color, 'MarkerSize', 15, 'LineWidth', leg_up_line_width);% leg 1 up
        Quadruped_leg1_down = plot(Quadruped_X(4:5, ii), Quadruped_Y(4:5, ii), '.-', 'Color', front_leg_color, 'MarkerSize', 15, 'LineWidth', leg_down_line_width);% leg 1 down
        hold on
        Quadruped_leg3_up = plot(Quadruped_X(6:7, ii), Quadruped_Y(6:7, ii), '.-', 'Color', front_leg_color, 'MarkerSize', 15, 'LineWidth', leg_up_line_width);% leg 3
        Quadruped_leg3_down = plot(Quadruped_X(7:8, ii), Quadruped_Y(7:8, ii), '.-', 'Color', front_leg_color, 'MarkerSize', 15, 'LineWidth', leg_down_line_width);% leg 3
        f = getframe(gcf);
        frames_video{ii} = f;
        if ii == 1
            [frames_gif(:,:,1,ii), cmap] = rgb2ind(f.cdata, 256, 'nodither');
        else
            frames_gif(:,:,1,ii) = rgb2ind(f.cdata, cmap, 'nodither');
        end
        if ii~=N_frames
            clf;
        end
    end
    % video file
    pathVideoAVI = [filename '.mp4']; % filename, used later to generate mp4
    writerObj = VideoWriter(pathVideoAVI,'MPEG-4');
    target_duration = (T/video_speed_up);
    writerObj.FrameRate = round(N_frames/target_duration);
    open(writerObj);
    for ii = 1:N_frames
        writeVideo(writerObj,frames_video{ii});
    end
    close(writerObj); % Close the movie file
end 
