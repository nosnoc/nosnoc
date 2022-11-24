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
%    License along with NOSNOC; if not, write to the Free Software Foundation,
%    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%

%% Info
% Simulation example from M. di Bernardo, K. H. Johansson, and F. Vasca. Self-oscillations and sliding in 
% relay feedback systems: Symmetry and bifurcations. International Journal of % Bifurcations and Chaos, 11(4):1121-1140, 2001
% and 
% Piiroinen, Petri T., and Yuri A. Kuznetsov. "An event-driven method to simulate Filippov systems with 
% accurate computing of sliding motions." ACM Transactions on Mathematical Software (TOMS) 34.3 (2008): 1-24.
%%
clear all
clc
close all
import casadi.*
%% discretization settings
N_finite_elements = 2;
T_sim = 10;
N_sim  = 200;

%% settings
settings = default_settings_nosnoc();
settings.use_fesd = 1;
settings.irk_scheme = 'Radau-IIA'; %'Gauss-Legendre';
settings.print_level = 2;
settings.n_s = 2;
settings.pss_mode = 'Stewart'; % 'Step;
settings.mpcc_mode = 3;  % Scholtes regularization
settings.comp_tol = 1e-9;
settings.cross_comp_mode  = 3;
settings.homotopy_update_rule = 'superlinear';
%% Time settings
model.N_finite_elements = N_finite_elements;
model.T_sim = T_sim;
model.N_sim = N_sim;
% Inital Value
model.x0 = [0;-0.001;-0.02];
% Variable defintion
x = SX.sym('x',3);
omega = 25;
xi = 0.05;
sigma = 1;

% xi = -0.07;
% sigma = 10;

A = [-(2*xi*omega+1)        1 0;...
     -(2*xi*omega+omega^2)  0 1;...
     -omega^2               0 0];
B = [1; -2*sigma;1];
c = x(1);
f_11 = A*x+B;
f_12 = A*x-B;

model.x = x;
model.c = c;
model.S = [-1;1];

F = [f_11 f_12];
model.F = F;
%% Call integrator
[results,stats,model] = integrator_fesd(model,settings);

%% Plot results
x1 = results.x_res(1,:);
x2 = results.x_res(2,:);
x3 = results.x_res(3,:);

if isequal(settings.pss_mode,'Stewart')
    theta = results.theta_res;
else
    alpha = results.alpha_res;
end
t_grid = results.t_grid;

figure
plot3(x2,x3,x1,'k-','LineWidth',1.5);
xlabel('$x_2$','Interpreter','latex');
ylabel('$x_3$','Interpreter','latex');
zlabel('$x_1$','Interpreter','latex');
grid on
% axis equal
% xlim([-4 4])
% ylim([-4 4])
% zlim([-4 4])
figure
subplot(121)
plot(t_grid,x1);
xlabel('$t$','Interpreter','latex');
ylabel('$x_1(t)$','Interpreter','latex');
grid on
subplot(122)
plot(t_grid,x1);
hold on
plot(t_grid,x2);
plot(t_grid,x3);
xlabel('$t$','Interpreter','latex');
ylabel('$x(t)$','Interpreter','latex');
grid on
legend({'$x_1(t)$','$x_2(t)$','$x_3(t)$'},'Interpreter','latex');

%%
figure
if isequal(settings.pss_mode,'Stewart')
    plot(t_grid,[[nan;nan],theta])
    xlabel('$t$','Interpreter','latex');
    ylabel('$\theta(t)$','Interpreter','latex');
    grid on    
    legend({'$\theta_1(t)$','$\theta_2(t)$'},'Interpreter','latex');
else
    plot(t_grid,[nan,alpha])
    xlabel('$t$','Interpreter','latex');
    ylabel('$\alpha(t)$','Interpreter','latex');
    grid on       
end