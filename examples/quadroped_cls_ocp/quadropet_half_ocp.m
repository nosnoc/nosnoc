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

%% Half quadroped
% example inspired by https://github.com/dojo-sim/ContactImplicitMPC.jl/tree/main/examples/quadruped
%%
clear all;
close all;
clc;
import casadi.*

% delete old gif
filename = 'half_quadroped_gait';
run_animation = 1;
video_speed_up = 1;
objective = 'jump'; % jump % wlk
%%
settings = NosnocOptions();
settings.irk_scheme = "RADAU_IIA";
settings.n_s = 2;  % number of stages in IRK methods
settings.use_fesd = 1;
settings.N_homotopy = 5;
settings.sigma_0 = 1e2;

settings.cross_comp_mode = 1;
settings.solver_opts.ipopt.tol = 1e-8;
settings.solver_opts.ipopt.acceptable_tol = 1e-8;
settings.solver_opts.ipopt.acceptable_iter = 3;

settings.solver_opts.ipopt.max_iter = 1000;

settings.equidistant_control_grid = 0;
settings.dcs_mode = 'CLS';

% settings.time_freezing = 1;
% settings.stabilizing_q_dynamics = 0;
% settings.s_sot_min = 0.99;
% settings.pss_lift_step_functions = 0;
% settings.stagewise_clock_constraint = 1;
% settings.nonsmooth_switching_fun = 0;

%% IF HLS solvers for Ipopt installed (check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions) use the settings below for better perfmonace:
settings.solver_opts.ipopt.linear_solver = 'ma57';
%% discretization
model = half_unitri_ai_model();
unfold_struct(model,'caller');

if isequal(objective,'walk')
    T = 2; % prediction horizon
    N_stg = 20; % control intervals
    N_FE = 3;  % integration steps per control intevral
    x_mid = x0;
    x_mid(1) = 0.75;
    x_end = x0;
    x_end(1) = 2;
    Q = diag([50; 50; 25*ones(n_q-2,1); 0.1*ones(n_q,1)])/10;
    Q_terminal= diag([500; 500; 500; 100*ones(n_q-3,1); 50*ones(n_q,1)])/10;
    u_ref = zeros(n_u,1);
    R = diag([0.1*ones(n_u,1)]);
elseif isequal(objective,'jump')
    T = 1; % prediction horizon
    N_stg = 10; % control intervals
    N_FE = 2;  % integration steps per control intevral
    x_mid = x0;
    x_mid(2) = 0.6;
    x_end = x0;
    x_end(2) = 0.8;
    Q = diag([10; 50;  10; 10; 10; 10; 10; 0.1*ones(n_q,1)]);
    Q_terminal = diag([300; 300; 300; 10; 10; 10; 10; 0.1*ones(n_q,1)]);
    Q_terminal = Q;
    u_ref = zeros(n_u,1);
    R = diag([0.01*ones(n_u,1)]);
    model.g_terminal = x(1:n_q)-x_end(1:n_q);
else
end
% interpolate refernece
x_ref = interp1([0 0.5 1],[x0,x_mid,x_end]',linspace(0,1,N_stg),'spline')'; %spline
%% Populate model
model.mu = 1;
model.T = T;
model.N_stages = N_stg;
model.N_finite_elements  = N_FE;

% model.u0 = 20*ones(n_u,1);
% model.a_n = 80;
% model.n_dim_contact = 2;

% LSQ objective
model.lsq_x = {x,x_ref,Q};
model.lsq_u = {u,u_ref,R};
model.lsq_T = {x,x_end,Q_terminal};

%% Call nosnoc solver
[results,stats,model,settings] = nosnoc_solver(model,settings);
%% read and plot results
unfold_struct(results,'base');
%
q_res = x_opt(1:n_q,:);
v_res = x_opt(n_q+1:2*n_q,:);

N_sim = length(q_res);
N_frames = N_sim;
h = T/(N_stg*N_FE);

if settings.time_freezing
    t_phy = x_opt(end,:);
    ind = diff(t_phy)/h>0.1;
    q_res = x_opt(1:n_q,ind);
    v_res = x_opt(n_q+1:2*n_q,ind);
    N_frames = sum(ind);
else
    N_frames = size(q_res,2);
end

if run_animation
    ground_X = [-0.1 2];
    ground_Y = [0  0];
    Quadruped_X = zeros(8, N_sim);
    Quadruped_Y = zeros(8, N_sim);
    for ii = 1:N_sim
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
    width = pos(3)*1.5;
    height = pos(4)*1.5;
    width = pos(3)*1.0;
    height = pos(4)*1.0;
    frames_gif = zeros(height, width, 1, N_frames, 'uint8');
    frames_video = { };
    for ii = 1:N_frames
        ground = plot(ground_X, ground_Y, '.-k', 'MarkerSize', 1, 'LineWidth', 2);
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