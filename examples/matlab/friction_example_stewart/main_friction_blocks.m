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

%% Info
% This is an example from 
% Stewart, D.E., 1996. A numerical method for friction problems with multiple contacts. The ANZIAM Journal, 37(3), pp.288-308.
% It considers 3 independent switching functions and it demonstrates the
% generalization of the FESD scheme presented in the NOS NOC software parep
%% settings
% collocation settings
[settings] = default_settings_fesd();  %% Optionally call this function to have an overview of all options.
settings.kappa = 0.05;
settings.n_s = 2;                            % Degree of interpolating polynomial
settings.irk_scheme = 'radau';     % Collocation scheme: radau or legendre
settings.mpcc_mode = 3;                    % 1 - extact, 2 - smooth  ,3 -relax , 4 - penalty, 5 - l1 penalty
settings.s_elastic_max = 1e1;                    
settings.cross_complementarity_mode = 8;
%% Generate Model
model = blocks_with_friction();
%% Simulation setings
model.T_sim = 12;
model.N_stages = 2;
model.N_finite_elements = 1;
model.h = 0.05;
settings.gamma_h = 1;
model.T = model.N_stages*model.h;
settings.use_previous_solution_as_initial_guess = 1;
%% Call FESD Integrator
[results,stats,model] = integrator_fesd(model,settings);
%% Formulate NLP; 
%% Get variables into main workspace
unfold_struct(model,'base');
unfold_struct(settings,'base');
unfold_struct(results,'base');
%% Read and plot Result 
% differential states
for i = 1:n_x
    eval( ['x' num2str(i) '_opt = results.x_res(' num2str(i) ',:);']);
end
for i = 1:n_theta
    eval( ['theta' num2str(i) '_opt = results.theta_res(' num2str(i) ',:);']);
end

for i = 1:n_theta
    eval( ['lambda' num2str(i) '_opt = results.lambda_res(' num2str(i) ',:);']);
end

for i = 1:n_simplex
    eval( ['mu' num2str(i) '_opt = results.mu_res(' num2str(i) ',:);']);
end
%%
figure
subplot(211)
plot(t_grid,x1_opt);
hold on
plot(t_grid,x2_opt);
plot(t_grid,x3_opt);
xlabel('$t$','interpreter','latex');
ylabel('$q(t)$','interpreter','latex');
legend({'$q_1(t)$','$q_2(t)$','$q_3(t)$'},'interpreter','latex');
grid on
subplot(212)
plot(t_grid,x4_opt);
hold on
plot(t_grid,x5_opt);
plot(t_grid,x6_opt);
xlabel('$t$','interpreter','latex');
ylabel('$v(t)$','interpreter','latex');
legend({'$v_1(t)$','$v_2(t)$','$v_3(t)$'},'interpreter','latex');
grid on
figure
plot(t_grid,x7_opt);
xlabel('$t$','interpreter','latex');
ylabel('$t(t)$','interpreter','latex');
grid on
% N_switches = sum((abs(diff(h_opt)))>1e-4);
% fprintf('Number of Switches is %d. \n',N_switches);
figure
subplot(311)
plot(t_grid,[theta1_opt,nan])
hold on
grid on
plot(t_grid,[theta2_opt,nan])
xlabel('$t$','interpreter','latex');
ylabel('$\theta(t)$','interpreter','latex');
legend({'$\theta_1(t)$','$\theta_2(t)$'},'interpreter','latex');
subplot(312)
plot(t_grid,[theta3_opt,nan])
hold on
grid on
plot(t_grid,[theta4_opt,nan])
xlabel('$t$','interpreter','latex');
ylabel('$\theta(t)$','interpreter','latex');
legend({'$\theta_3(t)$','$\theta_4(t)$'},'interpreter','latex');
subplot(313)
plot(t_grid,[theta5_opt,nan])
hold on
grid on
plot(t_grid,[theta6_opt,nan])
xlabel('$t$','interpreter','latex');
ylabel('$\theta(t)$','interpreter','latex');
legend({'$\theta_5(t)$','$\theta_6(t)$'},'interpreter','latex');

%%
figure
subplot(311)
plot(t_grid,[lambda1_opt,nan])
hold on
grid on
plot(t_grid,[lambda2_opt,nan])
xlabel('$t$','interpreter','latex');
ylabel('$\lambda(t)$','interpreter','latex');
legend({'$\lambda_1(t)$','$\lambda_2(t)$'},'interpreter','latex');
subplot(312)
plot(t_grid,[lambda3_opt,nan])
hold on
grid on
plot(t_grid,[lambda4_opt,nan])
xlabel('$t$','interpreter','latex');
ylabel('$\lambda(t)$','interpreter','latex');
legend({'$\lambda_3(t)$','$\lambda_4(t)$'},'interpreter','latex');
subplot(313)
plot(t_grid,[lambda5_opt,nan])
hold on
grid on
plot(t_grid,[lambda6_opt,nan])
xlabel('$t$','interpreter','latex');
ylabel('$\lambda(t)$','interpreter','latex');
legend({'$\lambda_5(t)$','$\lambda_6(t)$'},'interpreter','latex');
%%
figure
plot(t_grid,[mu1_opt,nan])
hold on
grid on
plot(t_grid,[mu2_opt,nan])
plot(t_grid,[mu3_opt,nan])
xlabel('$t$','interpreter','latex');
ylabel('$\mu(t)$','interpreter','latex');
legend({'$\mu_1(t)$','$\mu_2(t)$','$\mu_3(t)$'},'interpreter','latex');
%%
stairs(results.h_vec,'k')
xlabel('finite element','interpreter','latex');
ylabel('$h_{n}$','interpreter','latex');
%%
figure
semilogy(stats.complementarity_stats+1e-20,'k','LineWidth',1.5)
xlabel('integration step n','interpreter','latex');
ylabel('comp residual','interpreter','latex');
grid on

if model.N_stages == 2
    number_of_switches = sum(abs(diff(h_vec))>1e-6)
end