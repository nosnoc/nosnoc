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
clear; close all; clc;
import casadi.*

%%
play_animation = 1;

%% load nosnoc default options
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
%%
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 2;  % number of stages in IRK methods
problem_options.dcs_mode = 'CLS';
problem_options.cross_comp_mode = 7;
problem_options.friction_model = "Polyhedral";
problem_options.relax_terminal_numerical_time = 'ELL_2';
problem_options.rho_terminal_numerical_time = 1e3;
problem_options.print_level = 3;
problem_options.gamma_h = 0;
%problem_options.eps_cls = 0;

solver_options.sigma_0 = 1e0;
solver_options.homotopy_update_slope = 0.5;
solver_options.homotopy_update_rule = 'superlinear';
solver_options.N_homotopy = 100;
solver_options.print_level = 5;
% NLP solver settings;
default_tol = 1e-8;
solver_options.complementarity_tol = 1e-8;
solver_options.opts_casadi_nlp.ipopt.max_iter = 2e3;
solver_options.opts_casadi_nlp.ipopt.tol = default_tol;
solver_options.opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
solver_options.opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
solver_options.opts_casadi_nlp.ipopt.compl_inf_tol = default_tol;

%% IF HLS solvers for Ipopt installed use the settings below for better perfmonace (check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions) :
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';

%% discretizatioon
N_stg = 40; % control intervals
N_FE = 1;  % integration steps per control interval
T = 6;

problem_options.N_stages = N_stg;
problem_options.N_finite_elements  = N_FE;
problem_options.T = T;

%% model parameters
m1 = 1;
m2 = 1;
m3 = 1;
cart_width1 = 2;
cart_width2 = 2;
cart_width3 = 2;

M = diag([m1, m1, m2, m2, m3, m3]);

% Bounds on states and controls
ubx = ones(12,1)*10;
lbx = -ones(12,1)*10;
ubu = 30;
lbu = -30;

x0 = [-3; 1; 0; 1;  3; 1; ...
    0; 0; 0; 0; 0; 0];
u_ref = 0;

x_ref = [-7; 1; 0; 1; 5; 1;...
    0; 0; 0; 0; 0; 0];
 
Q = diag([10; 0.001; 1; 0.001; 10; 0.001; ...
          0.1; 0.1; 0.1; 0.1; 0.1; 0.1]);
Q_terminal = 100*Q;
R = 0.1;

%% Symbolic variables and bounds
g = 9.81;
q = SX.sym('q',6);
v = SX.sym('v',6);
u = SX.sym('u',1);
x = [q;v];

q1 = q(1:2);
q2 = q(3:4);
q3 = q(5:6);

model = nosnoc.model.Cls();
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

J_tangent = [0  0 1 0 0 ;...
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

D_tangent  = [];
for ii = 1:size(J_tangent,2)
    D_tangent  = [D_tangent, J_tangent(:,ii), -J_tangent(:,ii)];
end

model.f_c = f_c;
model.J_normal = J_normal;
model.J_tangent = J_tangent;
model.D_tangent = D_tangent;
model.e =  [0.0 0.5 0.0 0.0 0.0];
model.mu = [0.1 0.1 0.2 0.2 0.2];
% box constraints on controls and states
model.lbu = lbu;
model.ubu = ubu;
model.lbx = lbx;
model.ubx = ubx;
% Stage cost
model.f_q = (x-x_ref)'*Q*(x-x_ref)+ u'*R*u;
model.f_q_T = (x-x_ref)'*Q_terminal*(x-x_ref);

%% Call nosnoc solver
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
lambda_normal_guess = {};
ocp_solver.set('lambda_normal', 'init', {1:N_stg, 1:N_FE, 1:problem_options.n_s}, [0;0;g;g;g]);
ocp_solver.solve();

%% read and plot results
x = ocp_solver.get('x');
u_opt = ocp_solver.get('u');
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
t_grid = ocp_solver.get_time_grid();
t_grid_u = ocp_solver.get_control_grid();

%% animation
filename = 'three_carts_with_friction.gif';
x_min = min(x_ref)-2.5;
x_max = max(x_ref)+2.5;
cart_height = 2;
carts_appart = 1.5;

if play_animation
    figure(1)
    for ii = 1:length(p1x)
        % cart 1
        xp = [p1x(ii)-cart_width1/2 p1x(ii)+cart_height/2 p1x(ii)+cart_height/2 p1x(ii)-cart_width1/2];
        yp = [0 0 cart_height  cart_height];
        patch('XData',xp,'YData',yp,'FaceColor', [0 0.4470 0.7410],'FaceAlpha',0.9)
        hold on
        % cart 2
        xp = [p2x(ii)-cart_width2/2 p2x(ii)+cart_height/2 p2x(ii)+cart_height/2 p2x(ii)-cart_width2/2];
        yp = [0 0 cart_height  cart_height];
        patch('XData',xp,'YData',yp,'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha',0.9)

        % cart 3
        xp = [p3x(ii)-cart_width3/2 p3x(ii)+cart_height/2 p3x(ii)+cart_height/2 p3x(ii)-cart_width3/2];
        yp = [0 0 cart_height  cart_height];
        patch('XData',xp,'YData',yp,'FaceColor',[0.9290 0.6940 0.1250], 'FaceAlpha',0.9)
        %         % the refereneces
        % cart 1
        xp = [x_ref(1)-cart_width1/2 x_ref(1)+cart_height/2 x_ref(1)+cart_height/2 x_ref(1)-cart_width1/2];
        yp = [0 0 cart_height  cart_height];
        patch('XData',xp,'YData',yp,'FaceColor', [0 0.4470 0.7410],'FaceAlpha',0.1,'EdgeColor','none')
        hold on
        % cart 2
        xp = [x_ref(3)-cart_width2/2 x_ref(3)+cart_height/2 x_ref(3)+cart_height/2 x_ref(3)-cart_width2/2];
        yp = [0 0 cart_height  cart_height];
        patch('XData',xp,'YData',yp,'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha',0.1,'EdgeColor','none')

        % cart 3
        xp = [x_ref(5)-cart_width3/2 x_ref(5)+cart_height/2 x_ref(5)+cart_height/2 x_ref(5)-cart_width3/2];
        yp = [0 0 cart_height  cart_height];
        patch('XData',xp,'YData',yp,'FaceColor',[0.9290 0.6940 0.1250], 'FaceAlpha',0.1,'EdgeColor','none')

        % ground
        xp = [x_min x_max x_max x_min ];
        yp = [-1 -1 0 0];
        patch(xp,yp,0*ones(1,3),'FaceAlpha',0.1,'EdgeColor','none');

        axis equal
        xlim([x_min x_max])
        ylim([-0.75 3.5])
        pause(problem_options.h_k);

        frame = getframe(1);
        im = frame2im(frame);

        [imind,cm] = rgb2ind(im,256);
        if ii == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',problem_options.h_k(1));
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',problem_options.h_k(1));
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
    patch('XData',xp,'YData',yp,'FaceColor', [0 0.4470 0.7410],'FaceAlpha',0.9)
    hold on
    % cart 2
    xp = [p2x(ii)-cart_width2/2 p2x(ii)+cart_height/2 p2x(ii)+cart_height/2 p2x(ii)-cart_width2/2];
    yp = [0 0 cart_height  cart_height];
    patch('XData',xp,'YData',yp,'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha',0.9)

    % cart 3
    xp = [p3x(ii)-cart_width3/2 p3x(ii)+cart_height/2 p3x(ii)+cart_height/2 p3x(ii)-cart_width3/2];
    yp = [0 0 cart_height  cart_height];
    patch('XData',xp,'YData',yp,'FaceColor',[0.9290 0.6940 0.1250], 'FaceAlpha',0.9)

    %         % the refereneces
    % cart 1
    xp = [x_ref(1)-cart_width1/2 x_ref(1)+cart_height/2 x_ref(1)+cart_height/2 x_ref(1)-cart_width1/2];
    yp = [0 0 cart_height  cart_height];
    patch('XData',xp,'YData',yp,'FaceColor', [0 0.4470 0.7410],'FaceAlpha',0.2,'EdgeColor','none')
    hold on
    % cart 2
    xp = [x_ref(3)-cart_width2/2 x_ref(3)+cart_height/2 x_ref(3)+cart_height/2 x_ref(3)-cart_width2/2];
    yp = [0 0 cart_height  cart_height];
    patch('XData',xp,'YData',yp,'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha',0.2,'EdgeColor','none')

    % cart 3
    xp = [x_ref(5)-cart_width3/2 x_ref(5)+cart_height/2 x_ref(5)+cart_height/2 x_ref(5)-cart_width3/2];
    yp = [0 0 cart_height  cart_height];
    patch('XData',xp,'YData',yp,'FaceColor',[0.9290 0.6940 0.1250], 'FaceAlpha',0.2,'EdgeColor','none')

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

%%
lambda_tangent = ocp_solver.get("lambda_tangent");
lambda_tangent = [-lambda_tangent(5,:)+lambda_tangent(6,:);...
                    -lambda_tangent(7,:)+lambda_tangent(8,:);...
                  -lambda_tangent(9,:)+lambda_tangent(10,:);...
                   -lambda_tangent(1,:)+lambda_tangent(2,:);
                  -lambda_tangent(3,:)+lambda_tangent(4,:);
                 ];
lambda_normal = ocp_solver.get("lambda_normal");
lambda_normal = [lambda_normal(3:5,:);lambda_normal(1:2,:)];

figure('Renderer', 'painters', 'Position', [100 100 900 400])
% figure
subplot(231)
plot(t_grid,p1x,'LineWidth',1.5);
hold on
plot(t_grid,p2x,'LineWidth',1.5);
plot(t_grid,p3x,'LineWidth',1.5);
xlim([0 T])
% axis equal
grid on
legend({'$q_1(t)$','$q_2(t)$','$q_3(t)$'},'interpreter','latex','Location','southwest');
xlabel('$t$','interpreter','latex');
ylabel('$q(t)$','interpreter','latex');
% axis equal
subplot(232)
plot(t_grid,v1x,'LineWidth',1.5);
hold on
plot(t_grid,v2x,'LineWidth',1.5);
plot(t_grid,v3x,'LineWidth',1.5);
legend({'$v_1(t)$','$v_2(t)$','$v_3(t)$'},'interpreter','latex','Location','southeast');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v(t)$','interpreter','latex');
xlim([0 T])
subplot(233)
stairs(t_grid_u,[u_opt,nan],'LineWidth',1.5);
% legend({'$u_1(t)$','$u_2(t)$','$u_3(t)$'},'interpreter','latex');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$u(t)$','interpreter','latex');
xlim([0 T])

subplot(234)
Lambda_normal = ocp_solver.get("Lambda_normal");
stem(t_grid,[ones(2,1)*nan,Lambda_normal(1:2,:)]','LineWidth',1.5);
legend({'$\Lambda_{\mathrm{n}}^1(t)$','$\Lambda_{\mathrm{n}}^2(t)$'},'interpreter','latex','Location','northeast');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda_{\mathrm{n}}(t)$','interpreter','latex');
xlim([0 T])

subplot(235)
plot(t_grid,[ones(5,1)*nan,lambda_normal]','LineWidth',1.5);
legend({'$\lambda_{\mathrm{n}}^1(t)$','$\lambda_{\mathrm{n}}^2(t)$','$\lambda_{\mathrm{n}}^3(t)$','$\lambda_{\mathrm{n}}^4(t)$','$\lambda_{\mathrm{n}}^5(t)$'},'interpreter','latex','Location','east');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\lambda_{\mathrm{n}}(t)$','interpreter','latex');
ylim([-1 11])
xlim([0 T])

subplot(236)
plot(t_grid,[ones(size(lambda_tangent,1),1)*nan,lambda_tangent]','LineWidth',1.5);
legend({'$\lambda_{\mathrm{t}}^1(t)$','$\lambda_{\mathrm{t}}^2(t)$','$\lambda_{\mathrm{t}}^3(t)$','$\lambda_{\mathrm{t}}^4(t)$','$\lambda_{\mathrm{t}}^5(t)$'},'interpreter','latex','Location','southeast');
grid on
xlabel('$t$','interpreter','latex');
ylabel(['$\lambda_{\mathrm{t}}' ...
    '(t)$'],'interpreter','latex');

ylim([-2.5 2.5])
xlim([0 T])

