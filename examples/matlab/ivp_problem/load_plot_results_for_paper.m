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
%% Plot results for FESD paper

%%
all_scenarios_fesd = {'local_min_GL4_fesd_fixed','local_min_GL4_fesd_homotopy',...
                    'objective_GL4_fesd_fixed'};
all_scenarios_std =  {'local_min_GL4_std_fixed','local_min_GL4_std_homotopy',...
                    'objective_GL4_std_fixed'};



%% Objective function plots;

load([all_scenarios_std{3} '.mat'])
unfold_struct(results,'caller')
figure
plot(x0_vec,L_numeric,'LineWidth',2.0);
hold on
grid on
load([all_scenarios_fesd{3} '.mat'])
unfold_struct(results,'caller')
plot(x0_vec,L_numeric,'LineWidth',2.2);
plot(x0_vec,L_analytic,'k--','LineWidth',1.0);
xlabel('$x_0$','interpreter','latex');
ylabel('Objective','interpreter','latex');
legend({'$V_{\mathrm{Std}}(x_0)$','$V_{\mathrm{FESD}}(x_0)$','$V^{*}(x_0)$'},'interpreter','latex');

%% Error plots objective

figure
subplot(211)
semilogy(x0_vec,error_state)
grid on
ylabel('Error state - FESD','interpreter','latex');
subplot(212)
semilogy(x0_vec,error_objective)
grid on
xlabel('$x_0$','interpreter','latex');
ylabel('Error objective - FESD','interpreter','latex');

load([all_scenarios_std{3} '.mat'])
unfold_struct(results,'caller')

figure
subplot(211)
semilogy(x0_vec,error_state)
grid on
ylabel('Error state - Standard','interpreter','latex');
subplot(212)
semilogy(x0_vec,error_objective)
grid on
xlabel('$x_0$','interpreter','latex');
ylabel('Error objective - Standard','interpreter','latex');

%% Local minima expriments
matlab_blue = [0 0.4470 0.7410];
matlab_red = [0.8500 0.3250 0.0980];
load([all_scenarios_std{1} '.mat'])
unfold_struct(results,'caller')
figure
plot(x0_vec,x0_star,'color',matlab_blue,'LineWidth',2.0);
hold on
load([all_scenarios_std{2} '.mat'])
unfold_struct(results,'caller')
plot(x0_vec,x0_star,'color',matlab_blue,'LineStyle','-.','LineWidth',2.0);
load([all_scenarios_fesd{1} '.mat'])
unfold_struct(results,'caller')
plot(x0_vec,x0_star,'color',matlab_red,'LineWidth',2.0);
load([all_scenarios_fesd{2} '.mat'])
unfold_struct(results,'caller')
plot(x0_vec,x0_star,'color',matlab_red,'LineStyle','-.','LineWidth',2.0);
plot(x0_vec,x0_analytic,'k--','LineWidth',1.0);




grid on
axis equal
xlabel('$x_0$','interpreter','latex');
ylabel('$x_0^*$','interpreter','latex');
legend({'Analytic Solution','Standard (fixed)','Standard (Homotopy)','FESD (fixed)','FESD (Homotopy)'},'interpreter','latex');
