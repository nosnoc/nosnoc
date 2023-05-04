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

%
%
%% Quadroped OCP
% example inspired by https://github.com/dojo-sim/ContactImplicitMPC.jl/tree/main/examples/quadruped
%%
clear all;
close all;
clc;
import casadi.*

% delete old gif
filename = 'quadroped_jump_full';
run_animation = 1;
video_speed_up = 1;
objective = 'jump'; % jump % wlk
%%
settings = NosnocOptions();
settings.irk_scheme = 'RADAU_IIA';
settings.n_s = 2;  
settings.sigma_0 = 1e1;
settings.N_homotopy = 4;

settings.cross_comp_mode = 3;
settings.dcs_mode = 'CLS';

% settings.time_freezing = 1;
% settings.s_sot_min = 0.98;
% settings.equidistant_control_grid = 1;
% settings.pss_lift_step_functions = 0;
% settings.stagewise_clock_constraint = 1;
% settings.nonsmooth_switching_fun = 0;
% settings.stabilizing_q_dynamics = 1;

settings.solver_opts.ipopt.tol = 1e-8;
settings.solver_opts.ipopt.acceptable_tol = 1e-8;
settings.solver_opts.ipopt.acceptable_iter = 3;
settings.solver_opts.ipopt.max_iter = 700;

%% IF HLS solvers for Ipopt installed (check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions) use the settings below for better perfmonace:
settings.solver_opts.ipopt.linear_solver = 'ma27';

%% discretization
model = unitri_ai_model();
unfold_struct(model,'caller');


if isequal(objective,'walk')
    x0 = x0+0.05*rand(1);
    T = 1.5; % prediction horizon
    N_stg = 15; % control intervals
    N_FE = 2;  % integration steps per control intevral
    x_mid = x0;
    x_mid(1) = 1;
    x_mid(2) = 0.4;
%     x_mid(n_q+1) = 4;
    x_end = x0;
    x_end(1) = 1.5;
%     x_end(n_q+1) = 2;
    Q = diag([50; 10;  10; 10*ones(n_q-3,1); 10; 0.1*ones(n_q-1,1)]);
    Q = diag([1; 1;  1; 1*ones(n_q-3,1); 1*ones(n_q,1)]);
    Q_terminal= diag([500; 500; 500; 100*ones(n_q-3,1); 50*ones(n_q,1)]);
    u_ref = zeros(n_u,1);
    R = diag([0.1*ones(n_u,1)]);
elseif isequal(objective,'jump')
    T = 1.0; % prediction horizon
    N_stg = 10; % control intervals
    N_FE = 2;  % integration steps per control intevral
    x_mid = x0;
    x_mid(2) = 0.5;
    x_end = x0;
    x_end(2) = 0.7;
    Q = diag([100; 100;  100; 10*ones(n_q-3,1); 0.1*ones(n_q,1)]);
    Q_terminal = diag([300; 300; 300; 10*ones(n_q-3,1); 0.1*ones(n_q,1)]);
%     model.g_terminal = x(1:n_q)-x_end(1:n_q);
    u_ref = zeros(n_u,1);
    R = diag(0.1*ones(n_u,1));
end

%% Initial control guess
% control guess that cancels gravity
u0 =[-1.2859
   -2.2573
   -1.2859
   -2.2573
   -1.0505
   -4.4035
   -1.0505
   -4.4035];
model.u0 = u0;
%% interpolate refernece
x_ref = interp1([0 0.5 1],[x0,x_mid,x_end]',linspace(0,1,N_stg),'spline')'; %spline
%% Populate model
model.mu = 0.8;

model.T = T;
model.N_stages = N_stg;
model.N_finite_elements  = N_FE;

model.n_dim_contact = 2;
% LSQ objective
model.lsq_x = {x,x_ref,Q};
model.lsq_u = {u,u_ref,R};
model.lsq_T = {x,x_end,Q_terminal};

% model.a_n = 1e2;
% model.n_dim_contact = 2;
%% Call nosnoc solver
[results,stats,model,settings] = nosnoc_solver(model,settings);
%% read and plot results
unfold_struct(results,'base');
%%
q_res = x_opt(1:n_q,:);
v_res = x_opt(n_q+1:2*n_q,:);

h = T/(N_stg*N_FE);

if settings.time_freezing
    t_phy = x_opt(end,:);ind = diff(t_phy)/h>0.1;
    q_res = x_opt(1:n_q,ind);
    v_res = x_opt(n_q+1:2*n_q,ind);
    N_frames = sum(ind);
else
    N_frames  = size(q_res,2);
end


if run_animation
    ground_X = [-0.1 2];
    ground_Y = [0  0];
    Quadruped_X = zeros(14, N_frames);
    Quadruped_Y = zeros(14, N_frames);
    for ii = 1:N_frames
        if ii == 1
            x_input = x0;
        else
            x_input = x_opt(:, ii);
        end
        [Quadruped_X(:, ii), Quadruped_Y(:, ii)] = getQuadrupedConfiguration(model, x_input(1:11, 1));
    end
    %
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
        % back lags
        hold on
        Quadruped_leg2_up = plot(Quadruped_X(6:7, ii), Quadruped_Y(6:7, ii), '.-', 'Color', back_leg_color, 'MarkerSize', 15, 'LineWidth', leg_up_line_width);% leg 2
        Quadruped_leg2_down = plot(Quadruped_X(7:8, ii), Quadruped_Y(7:8, ii), '.-', 'Color', back_leg_color, 'MarkerSize', 15, 'LineWidth', leg_down_line_width);% leg 2
        hold on
        Quadruped_leg4_up = plot(Quadruped_X(12:13, ii), Quadruped_Y(12:13, ii), '.-', 'Color', back_leg_color, 'MarkerSize', 15, 'LineWidth', leg_up_line_width);% leg 4
        Quadruped_leg4_down = plot(Quadruped_X(13:14, ii), Quadruped_Y(13:14, ii), '.-', 'Color', back_leg_color, 'MarkerSize', 15, 'LineWidth', leg_down_line_width);% leg 4
        % torso
        hold on
        Quadruped_torso = plot(Quadruped_X(1:2, ii), Quadruped_Y(1:2, ii), '.-','Color', torso_color, 'MarkerSize', 15, 'LineWidth', torso_line_width);% torso
        % fron legs
        hold on
        Quadruped_leg1_up = plot(Quadruped_X(3:4, ii), Quadruped_Y(3:4, ii), '.-', 'Color', front_leg_color, 'MarkerSize', 15, 'LineWidth', leg_up_line_width);% leg 1 up
        Quadruped_leg1_down = plot(Quadruped_X(4:5, ii), Quadruped_Y(4:5, ii), '.-', 'Color', front_leg_color, 'MarkerSize', 15, 'LineWidth', leg_down_line_width);% leg 1 down
        hold on
        Quadruped_leg3_up = plot(Quadruped_X(9:10, ii), Quadruped_Y(9:10, ii), '.-', 'Color', front_leg_color, 'MarkerSize', 15, 'LineWidth', leg_up_line_width);% leg 3
        Quadruped_leg3_down = plot(Quadruped_X(10:11, ii), Quadruped_Y(10:11, ii), '.-', 'Color', front_leg_color, 'MarkerSize', 15, 'LineWidth', leg_down_line_width);% leg 3
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
    pathVideoAVI = [filename '.mp4']; 
    writerObj = VideoWriter(pathVideoAVI,'MPEG-4');
    target_duration = (T/video_speed_up);
    writerObj.FrameRate = round(N_frames/target_duration);
    open(writerObj);
    for ii = 1:N_frames
        writeVideo(writerObj,frames_video{ii});
    end
    % Close the movie file
    close(writerObj); 
end
