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

%% Hopper gait OCP
% example inspired by https://github.com/KY-Lin22/NIPOCPEC and https://github.com/thowell/motion_planning/blob/main/models/hopper.jl

%%
clear all;
close all;
clc;
import casadi.*


filename = 'hopper_long.gif';
delete hopper_long.gif

%%
[settings] = default_settings_nosnoc();
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 2;  % number of stages in IRK methods

settings.use_fesd = 1;

settings.N_homotopy = 6;
settings.cross_comp_mode = 1;
settings.opts_ipopt.ipopt.max_iter = 1e3;
settings.opts_ipopt.ipopt.tol = 1e-6;
settings.opts_ipopt.ipopt.acceptable_tol = 1e-6;
settings.opts_ipopt.ipopt.acceptable_iter = 3;

% settings.print_level = 2;
settings.comp_tol = 1e-9;
settings.time_freezing = 1;
% settings.s_sot_max = 2;
settings.s_sot_min = 1;
settings.equidistant_control_grid = 1;
settings.pss_lift_step_functions = 1;
settings.stagewise_clock_constraint = 1;
settings.g_ineq_at_fe = 1; % evaluate path constraint on every integration step
settings.g_ineq_at_stg = 1; % evaluate path constraint on every stage point

%% IF HLS solvers for Ipopt installed (check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions) use the settings below for better perfmonace:
% settings.opts_ipopt.ipopt.linear_solver = 'ma57';

% The methods and time-freezing refomulation are detailed in https://arxiv.org/abs/2111.06759
%% discretizatioon
T = 2; % prediction horizon
N_stg = 30; % control intervals
N_FE = 3;  % integration steps per control intevral

%% Hopper model
q = SX.sym('q', 4);
v = SX.sym('v', 4);
x = [q;v];
u = SX.sym('u', 3);
% state equations
mb = 1;  % body
ml = 0.1; % link
Ib = 0.25; % body
Il = 0.025; % link
mu = 1; % friction coeficient
g = 9.81;

% inertia matrix
M = diag([mb + ml, mb + ml, Ib + Il, ml]);
% coriolis and graviry
C = [0;(mb + ml)*g;0;0];
% Control input matrix
B = [0, -sin(q(3));...
    0, cos(q(3));...
    1, 0;...
    0, 1];

% constraint normal
f_c_normal = [0; 1; q(4)*sin(q(3)); -cos(q(3))];
% constraint tangent;
f_c_tangent = [1; 0; q(4)*cos(q(3)); sin(q(3))];

% tangential and normal velocity of the contact point
v_tangent = f_c_tangent'*v;
v_normal = f_c_normal'*v;

% All forces
f = -C + B*u(1:2);

% Gap function
f_c = q(2) - q(4)*cos(q(3));

% The control u(3) is a slack for modellingo of nonslipping constraints.
ubu= [50; 50; 100];
lbu= [-50; -50; 0];

g_comp_path = 0.01*[v_tangent*u(3);f_c*u(3)];

x0 = [0.1; 0.5; 0; 0.5; 0; 0; 0; 0];
ubx = [2; 1.5; pi; 0.5; 10; 10; 5; 5];
lbx =  [0; 0; -pi; 0.1; -10; -10; -5; -5];
        

Q = diag([50; 50; 20; 50; 0.1; 0.1; 0.1; 0.1]);
Q_terminal =diag([300; 300; 300; 300; 0.1; 0.1; 0.1; 0.1]);

u_ref = [0; 0; 0];
R = diag([0.001; 0.001;1e-5]);

%% interpolate refernece
    x_mid1 = [0.4; 0.65; 0; 0.2; 0; 0; 0; 0];
    x_mid2 = [0.6; 0.5; 0; 0.5; 0; 0; 0; 0];
    x_mid3 = [0.9; 0.65; 0; 0.2; 0; 0; 0; 0];
    x_end  = [1.3; 0.5; 0; 0.5; 0; 0; 0; 0];

    x_ref1 = interp1([0 0.5 1],[x0,x_mid1,x_mid2]',linspace(0,1,floor(N_stg/2)),'spline')'; %spline
    x_ref2 = interp1([0 0.5 1],[x_mid2,x_mid3,x_end]',linspace(0,1,ceil(N_stg/2)),'spline')'; %spline
    x_ref = [x_ref1,x_ref2];

%% Fill in model
model.T = T;
model.N_stages = N_stg;
model.N_finite_elements  = N_FE;
model.x = x;
model.u = u;
model.e = 0;
model.mu = mu;
model.a_n = 25;
% model.a_n = 1e2;
model.x0 = x0;

model.M = M;
model.f = f;
% gap functions
model.f_c = f_c;
model.tangent1 = f_c_tangent;
model.nabla_q_f_c = f_c_normal;

% box constraints on controls and states
model.lbu = lbu;
model.ubu = ubu;
model.lbx = lbx;
model.ubx = ubx;
% constraint on drift velocity
model.g_comp_path = g_comp_path;

model.lsq_x = {x,x_ref,Q};
model.lsq_u = {u,u_ref,R};
model.lsq_T = {x,x_end,Q_terminal};

%% Call nosnoc solver
[results,stats,model,settings] = nosnoc_solver(model,settings);

%% read and plot results
unfold_struct(results,'base');
q_opt = x_opt(1:4,:);
v_opt = x_opt(5:8,:);
t_opt = x_opt(9,:);
% gap function, normal and tangential velocity
c_eval = [];
for ii = 1:length(q_opt)
    % Note: Empty parameter argument as there are no global/time_varying params.
    c_eval = [c_eval,full(model.c_fun(x_opt(:,ii),[]))];
end

%%  plots
fig_num = 2;
plotHopperStatesControls(x_opt,u_opt,fig_num);
fig_num = 4;
% tagnetinal and normal velocity
plotHopperSwitchingFun(t_opt,c_eval,fig_num);
%%
if 0
figure
% slip complement
u_normalized = [u_opt(3,:),nan];
slip_cc1 = u_normalized .*c_eval(1,1:N_FE:end);
slip_cc2 = u_normalized .*c_eval(3,1:N_FE:end);
subplot(211)
plot(c_eval(1,1:N_FE:end));
hold on 
plot(u_normalized)
xlabel('$t$','interpreter','latex');
ylabel('Slip comelementarity','interpreter','latex');
legend({'Slack Gap','Gap'},'Interpreter','latex','Location','best')
subplot(212)
plot(slip_cc1)
hold on
plot(slip_cc2)
xlabel('$t$','interpreter','latex');
ylabel('Slip comelementarity','interpreter','latex');
legend({'Comp - slack vs gap','Comp - slack vs. v tan'},'Interpreter','latex','Location','best')
end
%% animation

x_min = x0(1)-0.4;
x_max = x_ref(1,end)+0.4;
y_min = -0.15;
y_max = 1;
q_ref = x_end(1:4);

x_head = q_opt(1,:);
y_head = q_opt(2,:);
x_foot = q_opt(1,:) - q_opt(4,:).*sin(q_opt(3,:));
y_foot = q_opt(2,:) - q_opt(4,:).*cos(q_opt(3,:));


figure(77)
for ii = 1:length(q_opt)
    % The reference
    plotHopperConfiguration(q_ref,0.1);
    hold on
    plot(x_ref(1,:),x_ref(2,:),'Color',[0 0 0 0.2]);
    % The Robot
    plotHopperConfiguration(q_opt(:,ii));
    % faded trajectory
    plot(x_head(1:ii),y_head(1:ii),'Color',[0 0 1 0.5]);
    hold on
    plot(x_foot(1:ii),y_foot(1:ii),'Color',[1 0 0 0.5]);
    plotHopperScene1(x_min,x_max,y_min,y_max);

    frame = getframe(77);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if ii == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',model.h_k(1));
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',model.h_k(1));
    end

    if ii~=length(q_opt)
        clf;
    end
end
