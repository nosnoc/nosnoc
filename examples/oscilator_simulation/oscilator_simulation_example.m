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

%%
clear all; clc; close all;
import casadi.*
%% basic settings
plot_integrator_output = 1;
plot_continious_time_sol = 1;
smooth_model = 0; % if 1, use model without switch for a sanity check
%% discretization settings
T_sim = pi/2;
N_sim  = 29;
N_FE = 2;

%% Init
problem_options = nosnoc.Options();
integrator_options = nosnoc.integrator.Options();
solver_options = integrator_options.fesd_solver_opts; % the fesd integrator uses an mpec solver, call and modify its options
model = nosnoc.model.Pss();


% select integrator
%integrator_options.integrator_plugin = "FESD";
integrator_options.integrator_plugin = "SMOOTHED_PSS";

% integraotr options
problem_options.use_fesd = 1;       % switch detection method on/off
problem_options.rk_scheme = RKSchemes.RADAU_IIA; %'Gauss-Legendre';
problem_options.n_s = 4;
problem_options.dcs_mode = 'Stewart'; % 'Step;
problem_options.N_finite_elements = N_FE; % number of finite elements
problem_options.T_sim = T_sim; % total simulation times 
problem_options.N_sim = N_sim; % number of simulations step during T_sim

% MPEC solver options
solver_options.print_level = 2;
solver_options.homotopy_steering_strategy = HomotopySteeringStrategy.ELL_INF;
solver_options.complementarity_tol = 1e-9; % Penalty/Relaxation parameter

%% Set up model;
x0 = [exp(-1);0]; % inital value
x_star = [exp(T_sim-1)*cos(2*pi*(T_sim-1));-exp((T_sim-1))*sin(2*pi*(T_sim-1))]; % Analytic solution, if T > 1;
% model parameters
omega = -2*pi;
R_osc  = 1;

A1 = [1 omega;...
    -omega 1];
A2 = [1 -omega;...
    omega 1];
if smooth_model
    A2 = A1;
end

% Variable defintion
x = SX.sym('x',2);
c = x'*x-1; % the switching surface is the unit circle
f_1 = A1*x;
f_2 = A2*x;
F = [f_1 f_2];

% Populate model
model.F = F;
model.x = x;
model.c = c;
model.S = [-1;1];
model.x0 = x0;

%% Call integrator
integrator = nosnoc.Integrator(model, problem_options, integrator_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();

%% numerical error
x_fesd = x_res(:,end);
error_x = norm(x_fesd-x_star,"inf");
fprintf(['Numerical error with h = %2.3f and ' char(problem_options.rk_scheme) ' with n_s = %d stages is: %5.2e: \n'],problem_options.h_sim,problem_options.n_s,error_x);
%% plot_solution_trajectory
t_star = R_osc; % eact switching time
if problem_options.use_fesd
    h_res = integrator.get("h");
    h_opt_full = h_res;
else
    h_opt_full = problem_options.h_k*ones(problem_options.N_sim*problem_options.N_stages,1);
end
x1_opt = x_res(1,:);
x2_opt = x_res(2,:);
if plot_integrator_output
    figure
    subplot(121)
    plot(t_grid,x1_opt,'linewidt',1.0);
    grid on
    hold on
    plot(t_grid,x2_opt,'linewidt',1.0);
    hh = -3:1:3;
    plot(hh*0+t_star,hh,'k')
    xlabel('$t$','interpreter','latex');
    ylabel('$x(t)$','interpreter','latex');
    legend({'$x_1(t)$','$x_2(t)$'},'interpreter','latex');
    subplot(122)
    plot(x1_opt,x2_opt,'linewidt',1.8);
    hold on
    theta = 0:0.01:2*pi;
    x = R_osc*(cos(theta));
    y = R_osc*(sin(theta));
    plot(x,y,'r','linewidth',1.5)
    xlabel('$x_1$','interpreter','latex');
    ylabel('$x_2$','interpreter','latex');
    axis equal
    grid on
    [X,Y] = meshgrid(-3:0.35:3,-3:0.35:3);
    [U,V] = oscilator_eval(X,Y);
    quiver(X,Y,U,V,'Color',0.65*ones(3,1),'LineWidth',1.2,'AutoScaleFactor',3);
    xlim([-exp(1) exp(1)]*1.01)
    ylim([-exp(1) exp(1)]*0.9)
    %
    figure
    stairs(h_opt_full)
    ylim([min(h_opt_full)*0.8 max(h_opt_full)*1.2])
    xlabel('integration step')
    ylabel('$h$','Interpreter','latex');
end

%% plot_continious_time_sol
if plot_continious_time_sol
    x_res_extended = x_res_full;
    tgrid_long = t_grid_full;

    %
    x1_very_fine = [];
    x2_very_fine = [];
    tgrid_very_fine = [];
    figure
    for ii =  1:problem_options.N_stages*problem_options.N_finite_elements*N_sim
        % read
        ind_now = 1+(ii-1)*(problem_options.n_s):(ii)*(problem_options.n_s)+1;
        tt = tgrid_long(ind_now);
        xx1 = x_res_extended(1,ind_now);
        xx2 = x_res_extended(2,ind_now);
        % fit
        p1 = polyfit(tt,xx1,length(xx1)-2);
        p2 = polyfit(tt,xx2,length(xx2)-2);
        t_eval = linspace(tt(1),tt(end),50);
        yy1 = polyval(p1,t_eval);
        yy2 = polyval(p2,t_eval);
        % store
        x1_very_fine = [x1_very_fine,yy1];
        x2_very_fine = [x2_very_fine,yy2];
        tgrid_very_fine = [tgrid_very_fine,t_eval];
        % plot
        plot(t_eval,yy1,'b');
        hold on;
        plot(t_eval,yy2,'r');
        plot(tt,xx1,'b.');
        plot(tt,xx2,'r.');
        grid on
    end
    xline(R_osc,'m')
    for ii = 1:length(t_grid)
        xline(t_grid(ii),'k--')
    end
end
