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
% Example 6.1 from 
% Dieci, Luca, and Luciano Lopez. 
% "Sliding motion in Filippov differential systems: theoretical results and a computational approach." 
% SIAM Journal on Numerical Analysis 47.3 (2009): 2023-2051.

%%
clear all
clc
close all
import casadi.*
%% discretization settings
N_finite_elements = 2;
T_sim = 40;
N_sim  = 100;

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
model.x0 = [0.04;-0.01];
model.x0 = [-0.8;-1];
% Variable defintion
x1 = SX.sym('x1',1);
x2 = SX.sym('x2',1);
x = [x1;x2];

c = x2-0.2;

f_11 = [x2;-x1+1/(1.2-x2)];
f_12 = [x2;-x1-1/(0.8+x2)];


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

if isequal(settings.pss_mode,'Stewart')
    theta = results.theta_res;
else
    alpha = results.alpha_res;
end
t_grid = results.t_grid;

figure
subplot(121)
plot(t_grid,x1);
xlabel('$t$','Interpreter','latex');
ylabel('$x_1(t)$','Interpreter','latex');
grid on
subplot(122)
plot(t_grid,x2);
xlabel('$t$','Interpreter','latex');
ylabel('$x_2(t)$','Interpreter','latex');
grid on

figure
plot(x1,x2)
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
grid on

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