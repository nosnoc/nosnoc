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
subplot(121)
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
% 
% figure
% subplot(211)
% semilogy(x0_vec,error_state)
% grid on
% ylabel('Error state - FESD','interpreter','latex');
% subplot(212)
% semilogy(x0_vec,error_objective)
% grid on
% xlabel('$x_0$','interpreter','latex');
% ylabel('Error objective - FESD','interpreter','latex');
% 
% load([all_scenarios_std{3} '.mat'])
% unfold_struct(results,'caller')
% 
% figure
% subplot(211)
% semilogy(x0_vec,error_state)
% grid on
% ylabel('Error state - Standard','interpreter','latex');
% subplot(212)
% semilogy(x0_vec,error_objective)
% grid on
% xlabel('$x_0$','interpreter','latex');
% ylabel('Error objective - Standard','interpreter','latex');

%% Local minima expriments
matlab_blue = [0 0.4470 0.7410];
matlab_red = [0.8500 0.3250 0.0980];
load([all_scenarios_std{1} '.mat'])
unfold_struct(results,'caller')
% figure
subplot(122)
plot(x0_vec,x0_star,'color',matlab_blue,'LineWidth',1.5);
hold on
load([all_scenarios_std{2} '.mat'])
unfold_struct(results,'caller')
plot(x0_vec,x0_star,'color',matlab_blue,'LineStyle','-.','LineWidth',1.5);
load([all_scenarios_fesd{1} '.mat'])
unfold_struct(results,'caller')
plot(x0_vec,x0_star,'color',matlab_red,'LineWidth',1.5);
load([all_scenarios_fesd{2} '.mat'])
unfold_struct(results,'caller')
plot(x0_vec,x0_star,'color',matlab_red,'LineStyle','-.','LineWidth',1.5);
plot(x0_vec,x0_analytic,'k-','LineWidth',1.0);
grid on
% axis equal
xlabel('$x_0$','interpreter','latex');
ylabel('$x_0^*$','interpreter','latex');
legend({'Standard (Fixed)','Standard (Homotopy)','FESD (Fixed)','FESD (Homotopy)','Analytic Solution'},'interpreter','latex');
