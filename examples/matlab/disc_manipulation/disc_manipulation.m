%
%    This file is part of NOSNOC.
%
%    NOSNOC -- A software for NOnSmooth Numerical Optimal Control.
%    Copyright (C) 2022 Armin Nurkanovic, Moritz Diehl (ALU Freiburg).
%
%    NOSNOC is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 3 of the License, or (at your option) any later version.
%
%    NOSNOC is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with NOS-NOC; if not, write to the Free Software Foundation,
%    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%

%% Manipulation of two discs

%%
clear all;
close all;
clc;
import casadi.*

%%
filename = 'discs_manipulation.gif';
%%
[settings] = default_settings_nosnoc();
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 1;  % number of stages in IRK methods

settings.use_fesd = 1;
settings.mpcc_mode = 5; % \ell_inifnity penalization of the complementariy constraints
settings.N_homotopy = 7;

settings.cross_comp_mode = 3;
settings.opts_ipopt.ipopt.max_iter = 1e3;
settings.print_level = 3;
settings.comp_tol = 1e-9;


settings.time_freezing = 1;


%% IF HLS solvers for Ipopt installed (check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions) use the settings below for better perfmonace:
% settings.opts_ipopt.ipopt.linear_solver = 'ma57';

%% discretizatioon
N_stg = 10; % control intervals
N_FE = 3;  % integration steps per control intevral
T = 2;

%% model parameters
m1 = 2;
m2 = 1;
r1 = 0.3;
r2 = 0.2;


ubx = [10; 10;10; 10; 5; 5; 5; 5];
lbx = -ubx;

q10 = [-1; -2];
q20 = [-1;-1];
v10 = [0;0];
v20 = [0;0];

q_target1 = [-1; 1];
q_target2 = [0; 0];

x0 = [q10;q20;v10;v20];

u_ref = [0;0];
x_ref = [q_target1;q_target2;zeros(4,1)];

Q = diag([5;5;10;10;0*ones(4,1)]);
R = diag([0.1 0.1]);
Q_terminal = 100*Q;

u_max = [20;20];
u_min = -u_max;

%% Symbolic variables and bounds
q = SX.sym('q',4);
v = SX.sym('v',4);
u = SX.sym('u',2);

q1 = q(1:2);
q2 = q(3:4);

x = [q;v];
model.T = T;
model.N_stages = N_stg;
model.N_finite_elements  = N_FE;
model.x = x;
model.u = u;
model.e = 0;
model.mu = 0.0;
model.a_n = 10;
model.x0 = x0;


model.M = diag([m1;m1;m2;m2]); % inertia/mass matrix;
model.f = [u;...
    zeros(2,1)];

% gap functions
model.c = [norm(q1-q2)^2-(r1+r2)^2];


% box constraints on controls and states
model.lbu = u_min;
model.ubu = u_max;
model.lbx = lbx;
model.ubx = ubx;
%% Stage cost
model.f_q = (x-x_ref)'*Q*(x-x_ref)+ u'*R*u;
model.f_q_T = (x-x_ref)'*Q_terminal*(x-x_ref);
%% Call nosnoc solver
[results,stats,model,settings] = nosnoc_solver(model,settings);
%% read and plot results
unfold_struct(results,'base');
p1 = x_opt(1,:);
p2 = x_opt(2,:);
p3 = x_opt(3,:);
p4 = x_opt(4,:);
v1 = x_opt(5,:);
v2 = x_opt(6,:);
v3 = x_opt(7,:);
v4 = x_opt(8,:);
t_opt = x_opt(9,:);

%% animation
figure('Renderer', 'painters', 'Position', [100 100 1000 800])

x_min =min([p1,p2,p3,p4])-1;
x_max = max([p1,p2,p3,p4])+1;

tt = linspace(0,2*pi,100);
x1 = r1*cos(tt);
y1 = r1*sin(tt);

x2 = r2*cos(tt);
y2 = r2*sin(tt);

for ii = 1:length(p1)
    plot(x1+p1(ii),y1+p2(ii),'k-','LineWidth',2);
    hold on
    plot(x2+p3(ii),y2+p4(ii),'r-','LineWidth',2);

    plot(x1+q_target1(1),y1+q_target1(2),'color',[0 0 0 0.6]);
    plot(x2+q_target2(1),y2+q_target2(2),'color',[1 0 0 0.6]);

    axis equal
    xlim([x_min x_max])
    ylim([x_min x_max])
    xlabel('$x$ [m]','Interpreter','latex');
    ylabel('$y$ [m]','Interpreter','latex');

    % save gif
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if ii == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',model.h_k(1));
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',model.h_k(1));
    end

    if ii~=length(p1)
        clf;
    end

end

%%
if 1
    figure('Renderer', 'painters', 'Position', [100 100 1400 600])
    subplot(311)
    plot(t_opt,v1,'LineWidth',1.5);
    hold on
    plot(t_opt,v2,'LineWidth',1.5);
    legend({'$v_1(t)$','$v_2(t)$'},'interpreter','latex');
    xlabel('$t$','interpreter','latex');
    ylabel('$v(t)$','interpreter','latex');
    grid on
    % axis equal
    subplot(312)
    plot(t_opt,v3,'LineWidth',1.5);
    hold on
    plot(t_opt,v4,'LineWidth',1.5);
    grid on
    legend({'$v_3(t)$','$v_4(t)$'},'interpreter','latex');
    xlabel('$t$','interpreter','latex');
    ylabel('$v(t)$','interpreter','latex');
    subplot(313)
    stairs(t_opt(1:N_FE:end),[u_opt,nan*ones(2,1)]','LineWidth',1.5);
    legend({'$u_1(t)$','$u_2(t)$'},'interpreter','latex');
    grid on
    xlabel('$t$','interpreter','latex');
    ylabel('$u$','interpreter','latex');
end