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

%% Three cart manipulation example
clear all;
close all;
clc;
import casadi.*

%%
play_animation = 1;
no_friction = 0;

%%
[settings] = NosnocOptions();
settings.irk_scheme = IRKSchemes.RADAU_IIA;
settings.n_s = 2;  % number of stages in IRK methods
settings.dcs_mode = 'CLS';
settings.cross_comp_mode = 3;
settings.friction_model = "Conic";
settings.conic_model_switch_handling = "Abs";
settings.print_level = 3;

if 1
    settings.sigma_0 = 1e0;
    settings.homotopy_update_slope = 0.2;
    settings.homotopy_update_rule = 'superlinear';
    settings.N_homotopy = 5;
    settings.opts_casadi_nlp.ipopt.max_iter = 2e3;
else
    settings.gamma_h = 0.9999;
    settings.sigma_0 = 1e0;
    settings.mpcc_mode = "Scholtes_ineq";
    settings.elastic_scholtes = 1;
    settings.homotopy_update_slope = 0.2;
    %settings.homotopy_update_rule = 'superlinear';
    settings.N_homotopy = 6;
    settings.opts_casadi_nlp.ipopt.max_iter = 5e3;
end

%% IF HLS solvers for Ipopt installed use the settings below for better perfmonace (check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions) :
settings.opts_casadi_nlp.ipopt.linear_solver = 'ma27';

%% discretizatioon
N_stg = 15; % control intervals
N_FE = 3;  % integration steps per control intevral
T = 5;

%% model parameters
m1 = 1;
m2 = 1;
m3 = 1;
cart_width1 = 2;
cart_width2 = 2;
cart_width3 = 2;

M = diag([m1, m1, m2, m2, m3, m3]);

ubx = [10;  inf; 10;    inf;  10;  inf;  10; 10; 10; 10; 10; 10];
lbx = [-10; -inf; -10; -inf; -10; -inf; -10; -10; -10; -10; -10;-10];


ubx = [10;  2;   10;    2;  10;  2;  10; 10; 10; 10; 10; 10];
lbx = [-10; -2; -10; -2; -10; -2; -10; -10; -10; -10; -10;-10];

x0 = [ -3; 1; 0; 1;  3; 1; ...
        0; 0; 0; 0; 0; 0];
x_ref = [-7; 0; 0; 0; 6; 0;...
          0; 0; 0; 0; 0; 0];
u_ref = 0;

Q = diag([10; 0.1; 1; 0.1; 10; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1]);
Q_terminal = 300*Q;
R = 0.1;

u_max = 30;
u_min = -30;

%% Symbolic variables and bounds
g = 9.81;
q = SX.sym('q',6);
v = SX.sym('v',6);
u = SX.sym('u',1);
x = [q;v];

q1 = q(1:2);
q2 = q(3:4);
q3 = q(5:6);

model = NosnocModel();
model.T = T;
model.dims.N_stages = N_stg;
model.dims.N_finite_elements  = N_FE;
model.x = x;
model.u = u;
model.x0 = x0;

model.M = M;
model.f_v = [ 0;...
    -m1*g;
    u;...
    -m2*g;
    0;...
    -m3*g];

% gap functions
f_c = [q2(1) - q1(1) - 0.5*cart_width2 - 0.5*cart_width1;...
    q3(1) - q2(1) - 0.5*cart_width3 - 0.5*cart_width2;...
    q1(2)-cart_width1/2;...
    q2(2)-cart_width2/2;...
    q3(2)-cart_width3/2;...
    ];

J_tangent = [ 0  0 1 0 0 ;...
             -1 -1 0 0 0;...
              0  0 0 1 0;...
              1  0 0 0 0;...
              0  0 0 0 1;...
              0  1 0 0 0];

J_tangent =   J_tangent./vecnorm(J_tangent);

J_normal = full(f_c.jacobian(q));
J_normal_fun = Function('J_normal_fun',{q},{J_normal});
J_normal = full(J_normal_fun(x0(1:6)))';
J_normal  = J_normal ./vecnorm(J_normal);


model.f_c = f_c;
model.J_normal = J_normal;
model.J_tangent = J_tangent;
model.D_tangent = [J_tangent, -J_tangent];

model.e =  [0.0 1.0 0.0 0.0 0.0];
if no_friction
    model.mu = 0.0;
else
%     model.mu = [0.0 0.0 0.001 0.001 0.001];
    model.mu = [0.2 0.2 0.2 0.2 0.2];
end
% box constraints on controls and states
model.lbu = u_min;
model.ubu = u_max;
model.lbx = lbx;
model.ubx = ubx;
% Stage cost
model.f_q = (x-x_ref)'*Q*(x-x_ref)+ u'*R*u;
model.f_q_T = (x-x_ref)'*Q_terminal*(x-x_ref);
% model.g_terminal = [x-x_ref];
%% Call nosnoc solver
solver = NosnocSolver(model, settings);
lambda_normal_guess = {};
y_gap_guess = {};
x_guess = {};
for ii = 1:N_stg
    lambda_normal_guess{ii} = [0;0;g;g;g];
%     lambda_tangent_guess{ii} = model.mu(1)*[0;0;g;g;g];
    y_gap_guess{ii} = [1;1;0;0;0];
    x_guess{ii} = x_ref;
end
solver.set('lambda_normal',lambda_normal_guess');
% solver.set('w0',rand(size(solver.problem.w)))
% solver.set('lambda_tangent',lambda_tangent_guess');

% use a previous solution to initialize
if isfile('results.mat')
    old_res = load('results.mat');
    names = fieldnames(old_res.extended);
    for k=1:numel(names)
        name = names{k};
        if name == 'x'
            val = reshape(old_res.extended.x(:,2:end), 1, []);
        else
            val = reshape(old_res.extended.(name), 1, []);
        end
        val = val';
        solver.set(name, val);
    end
end

%solver.set('y_gap',y_gap_guess');
%solver.set('Y_gap',y_gap_guess');
%solver.set('x', x_guess');
%solver.set('x_left_bp', x_guess');
[results,stats] = solver.solve();
%% read and plot results
unfold_struct(results,'base');
p1x = x(1,:);
p2x = x(3,:);
p3x = x(5,:);
p1y = x(2,:);
p2y = x(4,:);
p3y = x(6,:);
v1x = x(7,:);
v2x = x(9,:);
v3x = x(11,:);
v1y = x(8,:);
v2y = x(10,:);
v3y = x(12,:);

%% animation
% figure('Renderer', 'painters', 'Position', [100 100 1000 400])
filename = 'three_carts_with_friction.gif';
    carts_appart = 2;
    x_min = min(x_ref)-2.5;
    x_max = max(x_ref)+2.5;
    cart_height = 2;

    carts_appart = 1.5*1;
if play_animation
    figure(1)
    
    for ii = 1:length(p1x)
        % cart 1
        xp = [p1x(ii)-cart_width1/2 p1x(ii)+cart_height/2 p1x(ii)+cart_height/2 p1x(ii)-cart_width1/2];
        yp = [0 0 cart_height  cart_height];
        patch(xp,yp,'b','FaceAlpha',0.8)
        hold on
        % cart 2
        xp = [p2x(ii)-cart_width2/2 p2x(ii)+cart_height/2 p2x(ii)+cart_height/2 p2x(ii)-cart_width2/2];
        yp = [0 0 cart_height  cart_height];
        patch(xp,yp,'r','FaceAlpha',0.8)

        % cart 3
        xp = [p3x(ii)-cart_width3/2 p3x(ii)+cart_height/2 p3x(ii)+cart_height/2 p3x(ii)-cart_width3/2];
        yp = [0 0 cart_height  cart_height];
        patch(xp,yp,'m','FaceAlpha',0.8)

        % the refereneces
        % cart 1
        xp = [x_ref(1)-cart_width1/2 x_ref(1)+cart_height/2 x_ref(1)+cart_height/2 x_ref(1)-cart_width1/2];
        yp = [0 0 cart_height  cart_height];
        patch(xp,yp,'b','FaceAlpha',0.10)
        hold on
        % cart 2
        xp = [x_ref(3)-cart_width2/2 x_ref(3)+cart_height/2 x_ref(3)+cart_height/2 x_ref(3)-cart_width2/2];
        yp = [0 0 cart_height  cart_height];
        patch(xp,yp,'r','FaceAlpha',0.10)

        % cart 3
        xp = [x_ref(5)-cart_width3/2 x_ref(5)+cart_height/2 x_ref(5)+cart_height/2 x_ref(5)-cart_width3/2];
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

        if ii~=length(p1x)
            clf;
        end
    end
end
%%  several frames
figure('Renderer', 'painters', 'Position', [100 100 750 600])

x_min = min([x_ref])-2;
x_max = max([x_ref])+2;

N_total = length(p1x);
N_shots = 10;
N_skip = round(N_total/N_shots);
for jj= 1:N_shots
    if jj ~=N_shots
        ii = (jj-1)*(N_skip)+1;
    else
        ii = N_total;
    end
    subplot(N_shots/2,2,jj)
    % cart 1
    xp = [p1x(ii)-cart_width1/2 p1x(ii)+cart_height/2 p1x(ii)+cart_height/2 p1x(ii)-cart_width1/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'b','FaceAlpha',0.8)
    hold on
    % cart 2
    xp = [p2x(ii)-cart_width2/2 p2x(ii)+cart_height/2 p2x(ii)+cart_height/2 p2x(ii)-cart_width2/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'r','FaceAlpha',0.8)

    % cart 3
    xp = [p3x(ii)-cart_width3/2 p3x(ii)+cart_height/2 p3x(ii)+cart_height/2 p3x(ii)-cart_width3/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'m','FaceAlpha',0.8)

    % the refereneces
    % cart 1
    xp = [x_ref(1)-cart_width1/2 x_ref(1)+cart_height/2 x_ref(1)+cart_height/2 x_ref(1)-cart_width1/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'b','FaceAlpha',0.10)
    hold on
    % cart 2
    xp = [x_ref(3)-cart_width2/2 x_ref(3)+cart_height/2 x_ref(3)+cart_height/2 x_ref(3)-cart_width2/2];
    yp = [0 0 cart_height  cart_height];
    patch(xp,yp,'r','FaceAlpha',0.10)

    % cart 3
    xp = [x_ref(5)-cart_width3/2 x_ref(5)+cart_height/2 x_ref(5)+cart_height/2 x_ref(5)-cart_width3/2];
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
eval(['print -dpdf -painters ' ['manipulation_frames2'] ])

%%
figure('Renderer', 'painters', 'Position', [100 100 1100 260])
% figure
subplot(141)
plot(t_grid,p1x,'LineWidth',1.5);
hold on
plot(t_grid,p2x,'LineWidth',1.5);
plot(t_grid,p3x,'LineWidth',1.5);
xlim([0 T])
% axis equal
grid on
legend({'$q_1(t)$','$q_2(t)$','$q_3(t)$'},'interpreter','latex','Location','east');
xlabel('$t$','interpreter','latex');
ylabel('$q(t)$','interpreter','latex');
% axis equal
subplot(142)
plot(t_grid,v1x,'LineWidth',1.5);
hold on
plot(t_grid,v2x,'LineWidth',1.5);
plot(t_grid,v3x,'LineWidth',1.5);
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
stem(t_grid,[ones(2,1)*nan,Lambda_normal(1:2,:)]','LineWidth',1.5);
legend({'$\Lambda_{\mathrm{n}}^1(t)$','$\Lambda_{\mathrm{n}}^2(t)$'},'interpreter','latex','Location','northwest');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda_{\mathrm{n}}(t)$','interpreter','latex');
xlim([0 T])


set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
eval(['print -dpdf -painters ' ['carts_states'] ])
%%

figure
subplot(121)
plot(t_grid,[ones(5,1)*nan,lambda_normal]','LineWidth',1.5);
% legend({'$\Lambda_{\mathrm{n}}^1(t)$','$\Lambda_{\mathrm{n}}^2(t)$'},'interpreter','latex');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\lambda_{\mathrm{n}}(t)$','interpreter','latex');
subplot(122)
plot(t_grid,[ones(size(lambda_tangent,1),1)*nan,lambda_tangent]','LineWidth',1.5);
% legend({'$\Lambda_{\mathrm{n}}^1(t)$','$\Lambda_{\mathrm{n}}^2(t)$'},'interpreter','latex');
grid on
xlabel('$t$','interpreter','latex');
ylabel(['$\lambda_{\mathrm{t}}' ...
    '(t)$'],'interpreter','latex');
