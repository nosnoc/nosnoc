%
%    This file is part of NOS-NOC.
%
%    NOS-NOC -- A software for NOnSmooth Numerical Optimal Control.
%    Copyright (C) 2022 Armin Nurkanovic, Moritz Diehl (ALU Freiburg).
%
%    NOS-NOC is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 3 of the License, or (at your option) any later version.
%
%    NOS-NOC is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with NOS-NOC; if not, write to the Free Software Foundation,
%    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%
clear all
clc
close all
import casadi.*
%% NOS-NOC settings
[settings] = default_settings_nosnoc();  %% Optionally call this function to have an overview of all options.
% settings.n_s = 2;
% N_finite_elements = 6;

settings.n_s = 2;
N_finite_elements = 6;

settings.irk_representation = 'differential';
% settings.irk_scheme = 'Explicit-RK';
% settings.irk_scheme = 'Lobatto-IIIA';
% settings.irk_scheme = 'Gauss-Legendre';
% settings.irk_scheme = 'Lobatto-IIIC';
settings.irk_scheme = 'Radau-IIA';

settings.print_level = 3;
settings.use_fesd = 1;
settings.mpcc_mode = 5;
settings.comp_tol = 1e-12;
settings.equidistant_control_grid  = 1;

illustrate_regions  = 1;
terminal_constraint = 1;
linear_control = 1;

settings.step_equilibration = 1;
settings.rho_h = 1e1;
settings.heuristic_step_equilibration = 0;
settings.step_equilibration_mode = 3;
%% model equations
x_target = [-pi/6;-pi/4];

if linear_control
    v0  = [0;0];
else
    v0 = [];
end
model.x0 = [2*pi/3;pi/3;v0];
model.T = 4;
model.N_stages = 6;
model.N_finite_elements = N_finite_elements;
% Variable defintion
x1 = MX.sym('x1');
x2 = MX.sym('x2');

v1 = MX.sym('v1');
v2 = MX.sym('v2');
if linear_control
    x = [x1;x2;v1;v2];
else
    x = [x1;x2];
end
model.x = x;
% Control
u1 = MX.sym('u1');
u2 = MX.sym('u2');
model.u = [u1;u2];
if linear_control
    u_max = 10;
    model.lbu  = -u_max*ones(2,1);
    model.ubu  = u_max*ones(2,1);
else
    u_max = 2;
    model.lbu  = -u_max*ones(2,1);
    model.ubu  = u_max*ones(2,1);
end

% Switching Functions
p = 2; a = 0.15; a1 = 0;
b = -0.05; q = 3;

c1 = x1+a*(x2-a1)^p;
c2 = x2+b*x1^q;
S1 = [1;-1];
S2 = [1;-1];
model.c = {c1,c2};
model.S = {S1,S2};

%% Modes of the ODEs layers (for all  i = 1,...,n_simplex);
% part independet of the nonsmoothness
if linear_control
    f_11 = [-1+v1;0;u1;u2];
    f_12 = [1+v1;0;u1;u2];
    f_21 = [0;-1+v2;u1;u2];
    f_22 = [0;1+v2;u1;u2];
else
    f_11 = [-1+u1;0];
    f_12 = [1+u1;0];
    f_21 = [0;-1+u2];
    f_22 = [0;1+u2];
end
F1 = [f_11 f_12];
F2 = [f_21 f_22];
model.F = {F1,F2};

%% Objective
% model.f_q = u1^2+u2^2;
if linear_control
    model.f_q = 1*(v1^2+v2^2)+0*(u1^2+u2^2);
else
    model.f_q = u1^2+u2^2;
end
if terminal_constraint
    model.g_terminal = [x(1:2)-x_target(1:2)];
else
    model.f_q_T = 100*(x(1:2)-x_target(1:2))'*(x(1:2)-x_target(1:2));
end
%% Solve and plot
[results,stats,model,settings] = nosnoc_solver(model,settings);

u_opt = results.u_opt;
f_opt = full(results.f);

t_grid_optimizer = [results.t_grid];
x_res_optimizer = [results.x_opt];
%%
fprintf('Objective value %2.4f \n',f_opt);
f_star = 6.616653254750982;
fprintf('Error %2.4e \n',norm(f_opt - f_star));


x_res_integrator = [];
t_grid_integrator = [];
t_end = 0;

if 0
    model.T_sim = 4/6;
    model.N_sim = 12;
    model.N_stages = 1;
    model.N_finite_elements = 2;
    model.g_terminal = [];
    model.g_terminal_lb = [];
    model.g_terminal_ub = [];
    settings.mpcc_mode = 5;
    settings.irk_scheme = 'Radau-IIA';
    settings.n_s = 3;
    settings.use_fesd = 1;
    settings.print_level = 2;

    for ii =  1:6
        model.lbu = u_opt(:,ii);
        model.ubu = u_opt(:,ii);
        model.u0 = u_opt(:,ii);
        [results_integrator,stats,model] = integrator_fesd(model,settings);
        model.x0 = results_integrator.x_res(:,end);
        x_res_integrator = [x_res_integrator,results_integrator.x_res];
        t_grid_integrator = [t_grid_integrator, results_integrator.t_grid+t_end];
        t_end = t_grid_integrator(end);
    end

    x_end = x_res_integrator(:,end);
    if linear_control
        x_end = x_end(1:2);
    end
    fprintf('Error %2.4e \n',norm(x_target - x_end));
else
    tspan = [0 4/6];
    y0 = [2*pi/3;pi/3;0;0];
    sigma_int = settings.comp_tol;
    tol = sigma_int/10;
    options = odeset('RelTol', tol, 'AbsTol', tol/10);
    x_res_integrator = [];
    t_grid_integrator = [];
    for ii = 1:6
        [t,y_res] = ode15s(@(t,y) ...
            ((0.5*(1+tanh((y(1)+a*(y(2)-a1)^p)/sigma_int))))*([-1+y(3);0;u_opt(1,ii);u_opt(2,ii)])+(1-(0.5*(1+tanh((y(1)+a*(y(2)-a1)^p)/sigma_int))))*([1+y(3);0;u_opt(1,ii);u_opt(2,ii)])+...
            ((0.5*(1+tanh((y(2)+b*y(1)^q)/sigma_int))))*([0;-1+y(4);u_opt(1,ii);u_opt(2,ii)])+(1-(0.5*(1+tanh((y(2)+b*y(1)^q)/sigma_int))))*([0;1+y(4);u_opt(1,ii);u_opt(2,ii)])...
            ,tspan, y0,options);

        y0 = y_res(end,:)';
        x_res_integrator = [x_res_integrator,y_res'];
        t_grid_integrator = [t_grid_integrator,t'+(ii-1)*4/6];
    end
    x_end = x_res_integrator(1:2,end);
    fprintf('Error %2.4e \n',norm(x_target - x_end));
end

%%
if 0
figure
subplot(211)
plot(t_grid_optimizer,x_res_optimizer(1:2,:))
grid on
xlabel('$t$','interpreter','latex');
ylabel('$x(t)$ - optimizer','interpreter','latex');
subplot(212)
plot(t_grid_integrator,x_res_integrator(1:2,:))
grid on
xlabel('$t$','interpreter','latex');
ylabel('$x(t)$ - integrator','interpreter','latex');
%
figure
subplot(211)
plot(t_grid_optimizer,x_res_optimizer(3:4,:))
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v(t)$ - optimizer','interpreter','latex');
subplot(212)
plot(t_grid_integrator,x_res_integrator(3:4,:))
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v(t)$ - integrator','interpreter','latex');

%
figure
subplot(121)
plot(x_res_optimizer(1,:),x_res_optimizer(2,:),'LineWidth',2)
grid on
hold on
plot(x_target(1),x_target(2),'rx')
if illustrate_regions
    hold on
    t2 = -5:0.01:5;
    plot(-a*(t2-a1).^p,t2,'k')
    hold on
    t1 = t2;
    plot(t1,-b*t1.^q,'k')
    grid on
    axis equal
    xlim([-1.5 1.5])
    ylim([-1.5 1.5])
end
subplot(122)
plot(x_res_integrator(1,:),x_res_integrator(2,:),'LineWidth',2)
hold on
plot(x_target(1),x_target(2),'rx')
grid on
axis equal
if illustrate_regions
    hold on
%     p = 2; a = -0.1; a1 = 0;
%     b = -0.05; q = 3;
    t2 = -5:0.01:5;
    plot(-a*(t2-a1).^p,t2,'k')
    hold on
    t1 = t2;
    plot(t1,-b*t1.^q,'k')
    grid on
    axis equal
    xlim([-1.5 1.5])
    ylim([-1.5 1.5])
end
end
%% 
sliding_mode_plot_for_paper
