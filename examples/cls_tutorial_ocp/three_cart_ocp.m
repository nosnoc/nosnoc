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

%% Three cart manipulation example
% example taken from https://github.com/KY-Lin22/NIPOCPEC

%% 
clear all;
close all;
clc;
import casadi.*

% delete old gif
delete three_carts2.gif

%%
[settings] = NosnocOptions();  
settings.irk_scheme = IRKSchemes.GAUSS_LEGENDRE;
settings.n_s = 3;  % number of stages in IRK methods
settings.dcs_mode = 'CLS';

% settings.mpcc_mode = 'elastic_ineq'; % \ell_inifnity penalization of the complementariy constraints
settings.N_homotopy = 4;
settings.opts_casadi_nlp.ipopt.max_iter = 2e3;
settings.homotopy_update_slope = 0.2;
settings.cross_comp_mode = 1;

settings.homotopy_update_slope = 0.2;
settings.homotopy_update_rule = 'superlinear';
settings.N_homotopy = 6;
% settings.gamma_h = 0.999;
%% IF HLS solvers for Ipopt installed (check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions) use the settings below for better perfmonace:
settings.opts_casadi_nlp.ipopt.linear_solver = 'ma57';

%% discretizatioon
N_stg = 25; % control intervals
N_FE = 2;  % integration steps per control intevral
T = 5;

%% model parameters
m1 = 1;
m2 = 1;
m3 = 1;
cart_width1 = 2;
cart_width2 = 2;
cart_width3 = 2;
c_damping = 2;

ubx = [10; 15; 15; 5; 5; 5]; 
lbx = [-15; -15; -10; -5; -5; -5];            
        
x0 = [-3; 0; 3; 0; 0; 0];
x_ref = [-7; 0; 6; 0; 0; 0];
u_ref = 0;

Q = diag([10; 1; 10; 0.1; 0.1; 0.1]);
Q_terminal = 200*Q;
R = 0.1;

u_max = 30;
u_min = -30;

%% Symbolic variables and bounds
q = SX.sym('q',3);
v = SX.sym('v',3); 
u = SX.sym('u',1);
x = [q;v];
model.T = T;
model.N_stages = N_stg;
model.N_finite_elements  = N_FE;
model.x = x;
model.u = u;
model.e = [0 1];
model.mu = 0.0;
model.x0 = x0; 

model.M = diag([m1;m2;m3]); % inertia/mass matrix;
model.f_v = [-c_damping*v(1);...
           u-c_damping*v(2);...
           -c_damping*v(3)];

% gap functions
model.f_c = [q(2) - q(1) - 0.5*cart_width2 - 0.5*cart_width1;...
           q(3) - q(2) - 0.5*cart_width3 - 0.5*cart_width2];



% box constraints on controls and states
model.lbu = u_min;
model.ubu = u_max;
model.lbx = lbx;
model.ubx = ubx;
% Stage cost
model.f_q = (x-x_ref)'*Q*(x-x_ref)+ u'*R*u;
% terminal cost
model.f_q_T = (x-x_ref)'*Q_terminal*(x-x_ref);
% model.g_terminal = [x-x_ref];
%% Call nosnoc solver
solver = NosnocSolver(model, settings);
[results,stats] = solver.solve();
% [results] = polish_homotopy_solution(model,settings,results,stats.sigma_k); % (experimental, projects and fixes active set at solution and solves NLP)
%% read and plot results
unfold_struct(results,'base');
p1 = results.x(1,:);
p2 = results.x(2,:);
p3 = results.x(3,:);
v1 = results.x(4,:);
v2 = results.x(5,:);
v3 = results.x(6,:);
t_grid = results.t_grid;

%% animation
% figure('Renderer', 'painters', 'Position', [100 100 1000 400])
figure(1)
filename = 'three_carts2.gif';
carts_appart = 2;
x_min =min([p1,p2,p3])-2.5;
x_max = max([p1,p2,p3])+2.5;
cart_height = 2;

carts_appart = 1.5*1;
for ii = 1:length(p1)
    % cart 1
    xp = [p1(ii)-cart_width1/2 p1(ii)+cart_height/2 p1(ii)+cart_height/2 p1(ii)-cart_width1/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'b','FaceAlpha',0.8)
    hold on
    % cart 2
    xp = [p2(ii)-cart_width2/2 p2(ii)+cart_height/2 p2(ii)+cart_height/2 p2(ii)-cart_width2/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'r','FaceAlpha',0.8)

    % cart 3
    xp = [p3(ii)-cart_width3/2 p3(ii)+cart_height/2 p3(ii)+cart_height/2 p3(ii)-cart_width3/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'m','FaceAlpha',0.8)

    % the refereneces
    % cart 1
    xp = [x_ref(1)-cart_width1/2 x_ref(1)+cart_height/2 x_ref(1)+cart_height/2 x_ref(1)-cart_width1/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'b','FaceAlpha',0.10)
    hold on
    % cart 2
    xp = [x_ref(2)-cart_width2/2 x_ref(2)+cart_height/2 x_ref(2)+cart_height/2 x_ref(2)-cart_width2/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'r','FaceAlpha',0.10)

    % cart 3
    xp = [x_ref(3)-cart_width3/2 x_ref(3)+cart_height/2 x_ref(3)+cart_height/2 x_ref(3)-cart_width3/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'m','FaceAlpha',0.10)


    % ground     
    xp = [x_min x_max x_max x_min ];
    yp = [-1 -1 0 0];
    patch(xp,yp,0*ones(1,3),'FaceAlpha',0.1,'EdgeColor','none');
    
    axis equal
    xlim([x_min x_max])
    ylim([-0.75 3.5])
    pause(solver.model.h_k);
   
    frame = getframe(1);
    im = frame2im(frame);

    [imind,cm] = rgb2ind(im,256);
    if ii == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',solver.model.h_k(1));
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',solver.model.h_k(1));
    end

    if ii~=length(p1)
        clf;
    end
end

%%  several frames
figure('Renderer', 'painters', 'Position', [100 100 1000 600])

x_min = min([x(:)])-2;
x_max = max([x(:)])+2;

N_total = length(p1);
N_shots = 8;
N_skip = round(N_total/N_shots);
for jj= 1:N_shots
    if jj ~=N_shots
    ii = (jj-1)*(N_skip)+1;
    else
        ii = N_total;
    end
    subplot(N_shots/2,2,jj)
     % cart 1
    xp = [p1(ii)-cart_width1/2 p1(ii)+cart_height/2 p1(ii)+cart_height/2 p1(ii)-cart_width1/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'b','FaceAlpha',0.8)
    hold on
    % cart 2
    xp = [p2(ii)-cart_width2/2 p2(ii)+cart_height/2 p2(ii)+cart_height/2 p2(ii)-cart_width2/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'r','FaceAlpha',0.8)

    % cart 3
    xp = [p3(ii)-cart_width3/2 p3(ii)+cart_height/2 p3(ii)+cart_height/2 p3(ii)-cart_width3/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'m','FaceAlpha',0.8)

    % the refereneces
    % cart 1
    xp = [x_ref(1)-cart_width1/2 x_ref(1)+cart_height/2 x_ref(1)+cart_height/2 x_ref(1)-cart_width1/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'b','FaceAlpha',0.10)
    hold on
    % cart 2
    xp = [x_ref(2)-cart_width2/2 x_ref(2)+cart_height/2 x_ref(2)+cart_height/2 x_ref(2)-cart_width2/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'r','FaceAlpha',0.10)

    % cart 3
    xp = [x_ref(3)-cart_width3/2 x_ref(3)+cart_height/2 x_ref(3)+cart_height/2 x_ref(3)-cart_width3/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'m','FaceAlpha',0.10)

    text(-1.5,3,['$t = ' num2str(round(t_grid(ii),2)) '\ s$'],'interpreter','latex');
    xlabel('$x$ [m]','Interpreter','latex');
    if mod(jj,2)
        ylabel('$y$ [m]','Interpreter','latex');
    end
    % ground     
    xp = [x_min x_max x_max x_min ];
    yp = [-1 -1 0 0];
    patch(xp,yp,0*ones(1,3),'FaceAlpha',0.1,'EdgeColor','none');
    
    axis equal
    xlim([x_min x_max])
    ylim([-0.5 3.5])
end
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
eval(['print -dpdf -painters ' ['carts_frames'] ])


%%
figure('Renderer', 'painters', 'Position', [100 100 1100 260])
% figure
subplot(141)
plot(t_grid,p1,'LineWidth',1.5);
hold on
plot(t_grid,p2,'LineWidth',1.5);
plot(t_grid,p3,'LineWidth',1.5);
xlim([0 T])
% axis equal
grid on
legend({'$q_1(t)$','$q_2(t)$','$q_3(t)$'},'interpreter','latex','Location','east');
xlabel('$t$','interpreter','latex');
ylabel('$q(t)$','interpreter','latex');
% axis equal
subplot(142)
plot(t_grid,v1,'LineWidth',1.5);
hold on
plot(t_grid,v2,'LineWidth',1.5);
plot(t_grid,v3,'LineWidth',1.5);
legend({'$v_1(t)$','$v_2(t)$','$v_3(t)$'},'interpreter','latex','Location','south');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v(t)$','interpreter','latex');
xlim([0 T])
subplot(143)
stairs(t_grid(1:N_FE:end),[results.u,nan],'LineWidth',1.5);
% legend({'$u_1(t)$','$u_2(t)$','$u_3(t)$'},'interpreter','latex');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$u(t)$','interpreter','latex');
xlim([0 T])


subplot(144)
stem(t_grid,[ones(2,1)*nan,Lambda_normal]','LineWidth',1.5);
legend({'$\Lambda_{\mathrm{n}}^1(t)$','$\Lambda_{\mathrm{n}}^2(t)$'},'interpreter','latex','Location','northwest');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda_{\mathrm{n}}(t)$','interpreter','latex');
xlim([0 T])


set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
eval(['print -dpdf -painters ' ['carts_states'] ])