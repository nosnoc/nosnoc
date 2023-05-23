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

%% Hopper OCP
% example inspired by https://github.com/KY-Lin22/NIPOCPEC and https://github.com/thowell/motion_planning/blob/main/models/hopper.jl

% The methods and time-freezing refomulation are detailed in https://arxiv.org/abs/2111.06759
%%
clear all;
close all;
clc;
import casadi.*

% delete old gif
delete hopper_simple.gif

%%
[settings] = NosnocOptions();
settings.irk_scheme = IRKSchemes.RADAU_IIA;
settings.mpcc_mode = MpccMode.elastic_two_sided;
% settings.irk_representation = 'differential';
settings.n_s = 2;  % number of stages in IRK methods
settings.use_fesd = 1;
settings.N_homotopy = 7;
settings.homotopy_update_slope = 0.1;
settings.sigma_0 = 1e2;
settings.opts_casadi_nlp.ipopt.tol = 1e-6;
settings.opts_casadi_nlp.ipopt.acceptable_tol = 1e-6;
settings.opts_casadi_nlp.ipopt.acceptable_iter = 3;
settings.cross_comp_mode = 1;
settings.opts_casadi_nlp.ipopt.max_iter = 5e3;
settings.comp_tol = 1e-9;
settings.time_freezing = 0;
% settings.s_sot_max = 2;
% settings.s_sot_min = 1;
settings.equidistant_control_grid = 0;
%settings.use_speed_of_time_variables = 1;
%settings.local_speed_of_time_variable = 1;
settings.pss_lift_step_functions = 0;
settings.stagewise_clock_constraint = 0;
settings.dcs_mode = "CLS";
settings.friction_model = 'Conic';
settings.conic_model_switch_handling = 'Lp';
settings.g_path_at_fe = 0; % evaluate path constraint on every integration step
settings.g_path_at_stg = 0; % evaluate path constraint on every stage point
settings.nonsmooth_switching_fun = 0;
%settings.mpcc_mode = MpccMode.elastic_ineq;
settings.nlpsol = 'ipopt';
%settings.nlpsol = 'snopt';
%settings.opts_casadi_nlp.snopt.Major_feasibility_tolerance = 1e-3;
%settings.opts_casadi_nlp.snopt.Minor_feasibility_tolerance = 1e-3;

%% IF HLS solvers for Ipopt installed (check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions) use the settings below for better perfmonace:
% settings.opts_casadi_nlp.ipopt.linear_solver = 'ma57';

%% discretization
T = 1; % prediction horizon
N_stg = 20; % control intervals
N_FE = 3;  % integration steps per control intevral

v_slip_bound = 0.001;
full_comp = 1;
x_goal = 0.7;
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
mu = 0.45;
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
J_normal = [0; 1; q(4)*sin(q(3)); -cos(q(3))];
% constraint tangent;
J_tangent = [1; 0; q(4)*cos(q(3)); sin(q(3))];

% tangential and normal velocity of the contact point
v_tangent = J_tangent'*v;
v_normal = J_normal'*v;

% All forces
f_v = -C + B*u(1:2);

% Gap function
f_c = q(2) - q(4)*cos(q(3));

% The control u(3) is a slack for modellingo of nonslipping constraints.
ubu= [50; 50; 100];
lbu= [-50; -50; 0];

ubx = [x_goal+0.1; 1.5; pi; 0.50; 10; 10; 5; 5];
lbx =  [0; 0; -pi; 0.1; -10; -10; -5; -5];

x0 = [0.1; 0.5; 0; 0.5; 0; 0; 0; 0];
x_mid = [(x_goal-0.1)/2+0.1; 0.8; 0; 0.1; 0; 0; 0; 0];
x_end = [x_goal; 0.5; 0; 0.5; 0; 0; 0; 0];

Q = diag([50; 50; 20; 50; 0.1; 0.1; 0.1; 0.1]);
Q_terminal =diag([500; 500; 500; 500; 0.1; 0.1; 0.1; 0.1]);

Q = diag([50; 50; 20; 50; 0.1; 0.1; 0.1; 0.1]);
Q_terminal =diag([1000; 1000; 1000; 1000; 0.1; 0.1; 0.1; 0.1]);

u_ref = [0; 0; 0];
R = diag([0.01; 0.01; 1e-5]);

% Avoid slipping motion
g_comp_path = [v_tangent*u(3);f_c*u(3)];
g_comp_path = [v_tangent,u(3);f_c,u(3)];

%% interpolate refernece
x_ref = interp1([0 0.5 1],[x0,x_mid,x_end]',linspace(0,1,N_stg),'spline')'; %spline

%% Populate model
model = NosnocModel();
model.T = T;
model.N_stages = N_stg;
model.N_finite_elements  = N_FE;
model.x = x;
model.u = u;
model.e = 0;
model.mu = mu;
model.a_n = 25;
model.a_n = 1e2;
model.x0 = x0;

model.M = M;
model.f_v = f_v;
% gap functions
model.f_c = f_c;
model.J_tangent = J_tangent;

% box constraints on controls and states
model.lbu = lbu;
model.ubu = ubu;
model.lbx = lbx;
model.ubx = ubx;
% constraint on drift velocity
model.g_comp_path = g_comp_path;
% LSQ objective
model.lsq_x = {x,x_ref,Q};
model.lsq_u = {u,u_ref,R};
model.lsq_T = {x,x_end,Q_terminal};


%% create nosnoc solver
solver = NosnocSolver(model, settings);

%% Initialize solver with reference
x_init = num2cell(x_ref, 1)';
solver.set('x', x_init);
solver.set('x_left_bp', x_init);

%% Run solver
[results,stats] = solver.solve();
model = solver.model;
%% read and plot results
unfold_struct(results,'base');
q_opt = results.x(1:4,:);
v_opt = results.x(5:8,:);
% t_opt = x_opt(9,:);

if settings.time_freezing
% Swtiching functions: Gap function, normal velocity, tangentail velocity
c_eval = [];
for ii = 1:length(results.x)
    % Note: Empty parameter argument as there are no global/time_varying params.
    c_eval = [c_eval,full(model.c_fun(results.x(:,ii),[]))];
end
end
%%  plots
fig_num = 2;
plotHopperStatesControls(results.x,results.u,x_end,fig_num,settings,results);
fig_num = 4;
% tagnetinal and normal velocity
%plotHopperSwitchingFun(t_opt,c_eval,fig_num);

%%
if 0 
figure
% slip complement
u_normalized = [results.u(3,:),nan]./max(results.u(3,:));
slip_cc = u_normalized .*c_eval(1,1:N_FE:end);
plot(c_eval(1,1:N_FE:end));
hold on 
plot(u_normalized)
plot(slip_cc)
% yline(v_slip_bound,'r--')
% ylim([min(slip_cc)-0.01 max(min(slip_cc))+0.1])
xlabel('$t$','interpreter','latex');
ylabel('$v_{\mathrm{t}} \cdot (1-\mathrm{step}(f_c(q)))$ -tangential velocitiy','interpreter','latex');
end
%% animation
filename = 'hopper_simple.gif';

x_min = [x0(1)-0.4];
x_max = [x_ref(1,end)+0.4];
y_min = [-0.15];
y_max = [1];

x_head = q_opt(1,:);
y_head = q_opt(2,:);
x_foot = q_opt(1,:) - q_opt(4,:).*sin(q_opt(3,:));
y_foot = q_opt(2,:) - q_opt(4,:).*cos(q_opt(3,:));

x_head_ref = x_end(1);
y_head_ref = x_end(2);
x_foot_ref = x_end(1)- x_end(4).*sin(x_end(3));
y_foot_ref = x_end(2) - x_end(4).*cos(x_end(3));
q_ref = x_end(1:4);

figure(77)

for ii = 1:length(x_head)
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
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',solver.model.h_k(1));
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',solver.model.h_k(1));
    end

    if ii~=length(x_head)
        clf;
    end
end




