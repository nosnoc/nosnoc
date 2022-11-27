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
%%
clear all
clc
close all
import casadi.*
%%
plot_integrator_output = 1;
plot_continious_time_sol = 0;
%% discretization settings
T_sim = pi/2;
N_sim  = 29;
N_finite_elements = 2;
R_osc  = 1;

%% settings
% collocation settings
settings = default_settings_nosnoc();
settings.use_fesd = 1;       % switch detection method on/off
settings.irk_scheme = 'Radau-IIA'; %'Gauss-Legendre';
settings.print_level = 2;
settings.n_s = 4;
settings.pss_mode = 'Stewart'; % 'Step;
settings.mpcc_mode = 3;  % Scholtes regularization
% Penalty/Relaxation paraemetr
settings.comp_tol = 1e-9;
settings.cross_comp_mode  = 3;

%% Time settings
x_star = [exp(1);0];
T = T_sim;
x_star = [exp(T-1)*cos(2*pi*(T-1));-exp((T-1))*sin(2*pi*(T-1))];

model.N_finite_elements = N_finite_elements;
model.T_sim = T_sim;
model.N_sim = N_sim;
model.R_osc = R_osc;
smooth_model = 0; % if 1, use model without switch for a sanity check
omega = -2*pi;
A1 = [1 omega;...
    -omega 1];
A2 = [1 -omega;...
    omega 1];
if smooth_model
    A2 = A1;
end
% Inital Value
model.x0 = [exp(-1);0];
% Variable defintion
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1;x2];
c = x1^2+x2^2-1;
model.x = x;
model.c = c;
model.S = [-1;1];
f_11 = A1*x;
f_12 = A2*x;
F = [f_11 f_12];
model.F = F;
%% Call integrator
[results,stats,model] = integrator_fesd(model,settings);
%% numerical error
x_fesd = results.x_res(:,end);
error_x = norm(x_fesd-x_star,"inf");
fprintf(['Numerical error with h = %2.3f and ' settings.irk_scheme ' with n_s = %d stages is: %5.2e: \n'],model.h_sim,settings.n_s,error_x);
%% plot_solution_trajectory
t_star = R_osc; % eact switching time
x_res = results.x_res;
if isempty(results.h_vec)
    h_opt_full = h*ones(N_sim*N_stages,1);
else
    h_opt_full = results.h_vec;
end
x1_opt = x_res(1,:);
x2_opt = x_res(2,:);
t_grid = results.t_grid;
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
    unfold_struct(settings,'caller')
    x_res_extended = results.x_res_extended;
    h_opt = results.h_vec;
    [A_irk,b_irk,c_irk,order_irk] = generate_butcher_tableu(n_s,irk_scheme);
    t_grid = results.t_grid;
    tgrid_long = 0;
    h_grid_long = [];
    for ii  = 1:N_sim*N_stages*N_finite_elements;
        if use_fesd
            h_yet = h_opt(ii);
        else
            h_yet = model.h;
        end
        for jj = 1:n_s
            tgrid_long = [tgrid_long;t_grid(ii)+c_irk(jj)*h_yet];
        end
        tgrid_long = [tgrid_long;t_grid(ii)+h_yet];
    end

    %
    x1_very_fine = [];
    x2_very_fine = [];
    tgrid_very_fine = [];
    figure
    for ii =  1:N_stages*N_finite_elements*N_sim
        % read
        ind_now = 1+(ii-1)*(n_s+1):(ii)*(n_s+1)+1;
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
