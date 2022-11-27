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

%% Pendulum swing up with friction
% example taken from https://github.com/KY-Lin22/NIPOCPEC
%% Clear and import CasADi
clear all;
close all;
clc;
import casadi.*

% delete old gif
delete cart_pole_with_friction.gif
%% Build problem
import casadi.*
[settings] = default_settings_nosnoc();
% Choosing the Runge - Kutta Method and number of stages
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 2;
% MPCC Method
settings.mpcc_mode = 3;
settings.N_homotopy = 10;
settings.print_level = 3;

%% IF HLS solvers for Ipopt installed (check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions) use the settings below for better perfmonace:
% settings.opts_ipopt.ipopt.linear_solver = 'ma57';

%% Discretization parameters
model.N_stages = 50; % number of control intervals
model.N_finite_elements = 2; % number of finite element on every control intevral
model.T = 4;    % Time horizon

%% Model parameters and defintion
q = SX.sym('q', 2);
v = SX.sym('v', 2);
x = [q;v];
u = SX.sym('u', 1); % control

m1 = 1; % cart
m2 = 0.1; % link
link_length = 1;
g = 9.81;
% Inertia matrix
M = [m1 + m2, m2*link_length*cos(q(2));...
    m2 *link_length*cos(q(2)),  m2*link_length^2];
% Coriolis force
C = [0, -m2 * link_length*v(2)*sin(q(2));...
    0,   0];

% F_all (all forces) = Gravity+Control+Coriolis (+Friction)
f_all = [0;-m2*g*link_length*sin(x(2))]+[u;0]-C*v;

%  there is friction bewteen cart and ground
F_friction = 2;
% Dynamics with $ v > 0$
f_1 = [v;...
        inv(M)*(f_all-[F_friction;0])];
% Dynamics with $ v < 0$
f_2 = [v;...
        inv(M)*(f_all+[F_friction;0])];
F = [f_1,f_2];
% switching function (cart velocity)
c = v(1);

% specify initial and end state, cost ref and weight matrix
x0 = [1; 0/180*pi; 0; 0]; % start downwards
x_ref = [0; 180/180*pi; 0; 0]; % end upwards

Q = diag([1; 100; 1; 1]);
Q_terminal = diag([10; 100; 10; 20]);
Q_terminal = diag([100; 100; 10; 10]);
% Q_terminal = 10*Q;
R = 1;

% bounds
ubx = [5; 240/180*pi; 20; 20];
lbx = [-0.0; -240/180*pi; -20; -20];
u_max = 30;
u_ref = 0;

%% fill in model
model.F = F;
model.c = c;
model.lbx = lbx;
model.ubx = ubx;
model.x = x;
model.x0 =  x0;
model.u = u;
model.lbu = -u_max;
model.ubu = u_max;
% Sign matrix % f_1 for c=v>0, f_2 for c=v<0
model.S = [1;-1];

% Stage cost
if 1
    % directly via generic stage and terminal costs
    model.f_q = (x-x_ref)'*Q*(x-x_ref)+ (u-u_ref)'*R*(u-u_ref);
    % terminal cost
    model.f_q_T = (x-x_ref)'*Q_terminal*(x-x_ref);
    % model.g_terminal = x-x_target;
else
    % via least squares cost interace (makes time variable reference possible)
    model.lsq_x = {x,x_ref,Q};
    model.lsq_u = {u,u_ref,R};
    model.lsq_T = {x,x_ref,Q_terminal};
end
%% Solve OCP
[results,stats,model,settings] = nosnoc_solver(model,settings);
%% plots
% unfold structure to workspace of this script
unfold_struct(results,'base');
q1_opt = x_opt(1,:);
q2_opt= x_opt(2,:);
v1_opt= x_opt(3,:);
v2_opt= x_opt(4,:);


%% Animation
filename = 'cart_pole_with_friction.gif';
figure('Renderer', 'painters', 'Position', [100 100 1200 600])
cart_center = 0.125;
cart_width1 = 0.25;
cart_height = cart_center*2;
pole_X = [q1_opt',q1_opt'+(link_length)*cos(q2_opt'-pi/2)];
pole_Y = [cart_center+0*q1_opt',cart_center+link_length*sin(q2_opt'-pi/2)];
x_min =-3;
x_max = 3;
for ii = 1:length(q1_opt)
    % pole
    plot(pole_X(ii,:),pole_Y(ii,:),'k','LineWidth',3);
    hold on
    % tail
    plot(pole_X(1:ii,2),pole_Y(1:ii,2),'color',[1 0 0 0.5],'LineWidth',0.5);
    hold on
    % cart
    xp = [q1_opt(ii)-cart_width1/2 q1_opt(ii)+cart_height/2 q1_opt(ii)+cart_height/2 q1_opt(ii)-cart_width1/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'k','FaceAlpha',0.8)

    % targent
    % pole
    plot([x_ref(1),x_ref(1)+(link_length)*cos(x_ref(2)-pi/2)],[cart_center+0,cart_center+link_length*sin(x_ref(2)-pi/2)],'color',[0 0 0 0.1],'LineWidth',3);
    % cart
    xp = [x_ref(1)-cart_width1/2 x_ref(1)+cart_height/2 x_ref(1)+cart_height/2 x_ref(1)-cart_width1/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'k','FaceAlpha',0.1)

    % ground
    xp = [x_min x_max x_max x_min ];
    yp = [-2 -2 0 0];
    patch(xp,yp,0*ones(1,3),'FaceAlpha',0.1,'EdgeColor','none');

    axis equal
    xlim([x_min x_max])
    ylim([-1 2])
    text(-1.5,1.5,['Time: ' num2str(t_grid(ii),'%.2f') ' s'],'interpreter','latex','fontsize',15)
    %     try
    %         exportgraphics(gcf,'cart_pole_with_friction.gif','Append',true);
    %     catch
    %         disp('the simple gif function is avilable for MATLAB2022a and newer.')
    %     end

    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if ii == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',model.h_k(1));
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',model.h_k(1));
    end

    if ii~=length(q1_opt)
        clf;
    end
end

%% states

figure
subplot(131)
plot(t_grid,q1_opt)
hold on
plot(t_grid,q2_opt)
yline(pi,'k--')
ylabel('$q(t)$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
grid on
legend({'$q_1(t)$ - cart ','$q_2(t)$ - pole'},'Interpreter','latex','Location','best')
subplot(132)
plot(t_grid,v1_opt)
hold on
plot(t_grid,v2_opt)
yline(0,'k--')
ylabel('$v(t)$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
grid on
legend({'$v_1(t)$ - cart ','$v_2(t)$ - pole'},'Interpreter','latex','Location','best')
% t_grid_u = t_grid_u';
subplot(133)
stairs(t_grid_u,[u_opt,nan]);
ylabel('$u(t)$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
grid on
xlim([0 model.T])
