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
function [varargout] = plot_result_ball(model,settings,results,stats)

unfold_struct(model,'caller');
unfold_struct(settings,'caller')
d = 3;

obj = full(results.f);
w_opt = full(results.x);

%  Colors
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
organe = [0.9290 0.6940 0.1250];
grey = [0.85 0.85 0.85];
%%
n = n_x+n_u;
nn = n-1;
tgrid = linspace(0, T, N_stages+1);
tgrid_z = linspace(0, T, N_stages);
%% read resultsutions

diff_states = w_opt(ind_x);
controls = w_opt(ind_u);
alg_states = w_opt(ind_z);

% differential states
for i = 1:n_x
    eval( ['x' num2str(i) '_opt = diff_states(' num2str(i) ':n_x+n_x*d:end);']);
end


% convex multiplers
for i = 1:n_theta
    eval( ['theta' num2str(i) '_opt = alg_states(' num2str(i) ':n_z+n_z*(d-1):end);']);
    %     eval( ['theta' num2str(i) '_opt = alg_states(' num2str(i) ':n_z:end);']);
end
% lambdas
for i = 1:n_theta
    %     eval( ['lambda' num2str(i) '_opt = alg_states(' num2str(i+n_theta) ':n_z+n_z*d:end);']);
    eval( ['lambda' num2str(i) '_opt = alg_states(' num2str(i+n_theta) ':n_z+n_z*(d-1):end);']);
end
% mu
for i = 1:n_simplex
    %     eval( ['mu_opt = alg_states(' num2str(i+2*n_theta) ':n_z+n_z*d:end);']);
    eval( ['mu' num2str(i) '_opt = alg_states(' num2str(i+2*n_theta) ':n_z+n_z*(d-1):end);']);
end

for i = 1:n_u
    %     eval( ['mu_opt = alg_states(' num2str(i+2*n_theta) ':n_z+n_z*d:end);']);
    eval( ['u' num2str(i) '_opt = controls(' num2str(i) ':n_u:end);']);
end

if use_fesd
    h_opt = w_opt(ind_h);
    tgrid = (cumsum([0;h_opt]));
    tgrid_z = cumsum(h_opt)';
end


%%
if mpcc_mode == 4
    ind_t = find([1;theta1_opt]>1e-2);
else
    ind_t = find(diff([nan;x5_opt;nan])>1e-5);
end

if 0
time_physical = x5_opt(ind_t);
% Geomtric plot
x_target = 10;
figure;
x = [0 x_target x_target 0];
y = [0 0 -1 -1];
patch(x,y,'k','FaceAlpha',0.2)
hold on
plot(x1_opt(ind_t),x2_opt(ind_t),'linewidth',1.2,'color',0*ones(3,1));
grid on
hold on

xlabel('$q_1$','interpreter','latex');
ylabel('$q_2$','interpreter','latex');
axis equal
ylim([-0.4 max(x2_opt)*1.15])
xlim([0.0 x_target])
saveas(gcf,'geometric_traj')
%

matlab_blue = [0 0.4470 0.7410];
matlab_red = [0.8500 0.3250 0.0980];
figure
subplot(121)
plot(x5_opt,x3_opt,'linewidth',1.2,'color',matlab_blue,'LineStyle','--');
hold on
plot(x5_opt,x4_opt,'linewidth',1.2,'color',matlab_red);

xlabel('$t$','interpreter','latex');
ylabel('$v(t)$','interpreter','latex');
grid on
legend({'$v_1(t)$','$v_2(t)$'},'interpreter','latex');
xlim([0 T]);
subplot(122)
stairs(x5_opt(1:N_finite_elements:end),[u1_opt;nan],'color',matlab_blue,'linewidth',1.2,'LineStyle','--');
hold on
stairs(x5_opt(1:N_finite_elements:end),[u2_opt;nan],'color',matlab_red,'linewidth',1.2);
grid on
xlabel('$t$','interpreter','latex');
ylabel('$u(t)$','interpreter','latex');
grid on
legend({'$u_1(t)$','$u_2(t)$'},'interpreter','latex');
xlim([0 T]);
%
saveas(gcf,'velocity_and_control')
end
%%
% Geomtric plot
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultlegendinterpreter','latex')
set(groot,'DefaultTextarrowshapeInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');
matlab_blue = [0 0.4470 0.7410];
matlab_red = [0.8500 0.3250 0.0980];
x_target = 4;
figure;
subplot(131)
x = [0 x_target x_target 0];
y = [0 0 -1 -1];
patch(x,y,'k','FaceAlpha',0.2)
hold on
plot(x1_opt(ind_t),x2_opt(ind_t),'linewidth',1.2,'color',0*ones(3,1));
grid on
hold on
xlabel('$q_1$','interpreter','latex');
ylabel('$q_2$','interpreter','latex');
% axis equal
ylim([-0.4 max(x2_opt)*1.15])
xlim([0.0 x_target])
subplot(132)
plot(x5_opt,x3_opt,'linewidth',1.2,'color',matlab_blue);
hold on
plot(x5_opt,x4_opt,'linewidth',1.2,'color',matlab_red);
xlabel('$t$','interpreter','latex');
ylabel('$v(t)$','interpreter','latex');
grid on
legend({'$v_1(t)$','$v_2(t)$'},'interpreter','latex','NumColumns',2);
xlim([0 T]);
ylim([-6 6])
subplot(133)
stairs(x5_opt(1:N_finite_elements:end),[u1_opt;nan],'color',matlab_blue,'linewidth',1.2);
hold on
stairs(x5_opt(1:N_finite_elements:end),[u2_opt;nan],'color',matlab_red,'linewidth',1.2);
grid on
xlabel('$t$','interpreter','latex');
ylabel('$u(t)$','interpreter','latex');
grid on
legend({'$u_1(t)$','$u_2(t)$'},'interpreter','latex','NumColumns',2);
xlim([0 T]);
ylim([-10 10])
%%

%% spee of time plots
if 0
    if use_fesd
        h_opt_stagewise = reshape(h_opt,N_finite_elements(1),N_stages);
        if length(ind_tf)>1
            s_sot = w_opt(ind_tf);
        elseif length(ind_tf) == 0
            s_sot = ones(N_stages,1);
        else
            s_sot = w_opt(ind_tf)*ones(N_stages,1);
        end
        t_control_grid_pseudo = cumsum([0,sum(h_opt_stagewise)]);
        t_control_grid_pseudo_streched = cumsum([0,sum(h_opt_stagewise).*s_sot']);
    end
    %%
    if use_fesd
        figure
        stairs(tgrid(1:N_finite_elements:end),[s_sot;nan],'linewidth',1.2)
        xlabel('$\tau$','interpreter','latex');
        ylabel('$s(\tau)$','interpreter','latex');
        grid on
        ylim([0.5 3.2])
    end

    figure
    tt = linspace(0,T,N_stages*N_finite_elements(1)+1);
    plot(tt,x5_opt,'linewidth',1.2)
    hold on
    plot(tt,x5_opt,'.','color',blue,'MarkerSize',7)
    plot(tt(1:N_finite_elements:end),x5_opt(1:N_finite_elements:end),'k.','MarkerSize',9)
    plot(tt,tt,'k:')
    % grid on
    h = T/N_stages;
    for ii = 0:N_stages+1
        xline(h*ii,'k--')
    end

    axis equal
    xlim([0 T])
    ylim([0 T])
    xlabel('$\tau$ [Numerical Time]','interpreter','latex');
    ylabel('$t(\tau)$ [Phyisical Time]','interpreter','latex');
end
end

