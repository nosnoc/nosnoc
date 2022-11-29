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
% Example 4.1 from 
% Dieci, Luca, and Nicola Guglielmi. "Regularizing piecewise smooth differential systems: co-dimension $ $2 $$ discontinuity surface." Journal of Dynamics and Differential Equations 25.1 (2013): 71-94.

%%
clear all
clc
close all
import casadi.*
%% discretization settings
N_finite_elements = 3;
T_sim = 20;
N_sim  = 100;

%% settings
settings = default_settings_nosnoc();
settings.use_fesd = 1;
settings.irk_scheme = 'Radau-IIA'; %'Gauss-Legendre';
settings.print_level = 2;
settings.n_s = 4;
settings.pss_mode = 'Stewart'; % 'Step;
settings.mpcc_mode = 3;  % Scholtes regularization
settings.comp_tol = 1e-9;
settings.cross_comp_mode  = 3;
settings.homotopy_update_rule = 'superlinear';
settings.pss_lift_step_functions = 0;
%% Time settings
model.N_finite_elements = N_finite_elements;
model.T_sim = T_sim;
model.N_sim = N_sim;
% Inital Value
model.x0 = [0.04;-0.01;-0.02];
% Variable defintion
x1 = SX.sym('x1',1);
x2 = SX.sym('x2',1);
x3 = SX.sym('x3',1);
x = [x1;x2;x3];

c = [x2-0.2; x3-0.4];

f_11 = [(x2+x3)/2;-x1+1/(1.2-x2);-x1+1/(1.4-x3)];
f_12 = [(x2+x3)/2;-x1+1/(1.2-x2);-x1+1/(0.6+x3)];
f_13 = [(x2+x3)/2;-x1-1/(0.8+x2);-x1+1/(1.4-x3)];
f_14 = [(x2+x3)/2+x1*(x2+0.8)*(x3+0.6);-x1-1/(0.8+x2);-x1-1/(0.6+x3)];

model.x = x;
model.c = c;
model.S = [-1 -1;-1 1;1 -1; 1 1];

F = [f_11 f_12 f_13 f_14];
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
subplot(131)
plot(t_grid,x1);
xlabel('$t$','Interpreter','latex');
ylabel('$x_1(t)$','Interpreter','latex');
grid on
subplot(132)
plot(t_grid,x2);
xlabel('$t$','Interpreter','latex');
ylabel('$x_2(t)$','Interpreter','latex');
grid on
subplot(133)
plot(t_grid,x3);
xlabel('$t$','Interpreter','latex');
ylabel('$x_3(t)$','Interpreter','latex');
grid on

figure
plot(x1,x3)
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_3$','Interpreter','latex');
grid on

%%
figure
if isequal(settings.pss_mode,'Stewart')
    plot(t_grid,[[nan;nan;nan;nan],theta])
    xlabel('$t$','Interpreter','latex');
    ylabel('$\theta(t)$','Interpreter','latex');
    grid on    
    legend({'$\theta_1(t)$','$\theta_2(t)$','$\theta_3(t)$','$\theta_4(t)$'},'Interpreter','latex');
else
    plot(t_grid,[[nan;nan],alpha])
    xlabel('$t$','Interpreter','latex');
    ylabel('$\alpha(t)$','Interpreter','latex');
    grid on       
end