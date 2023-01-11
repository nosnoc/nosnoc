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

x_opt_s = cellfun(@(x) w_opt(x), structured_ind.x, 'uni', 0);
theta_opt_s = cellfun(@(theta) w_opt(theta), structured_ind.theta, 'uni', 0);
lambda_opt_s = cellfun(@(lam) w_opt(lam), structured_ind.lam, 'uni', 0);
mu_opt_s = cellfun(@(mu) w_opt(mu), structured_ind.mu, 'uni', 0);

% differential states
x_opt = cell(n_x, 1);
x_opt_flat = cell(n_x, 1);
for i = 1:n_x
    x_opt{i} = cellfun(@(x) x(i), x_opt_s);
    x_opt_flat{i} = reshape(x_opt{i}, prod(size(x_opt{i})), 1);
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
for i = 1:n_sys
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
if strcmp(mpcc_mode, 'Sholtes_eq')
    ind_t = find([1;theta1_opt]>1e-2);
else
    ind_t = find(diff([nan;x_opt_flat{1};nan])>1e-5);
end

% TODO: fix this
if 0
    time_physical = x5_opt(ind_t);
    % Geomtric plot
    x_target = 10;
    figure;
    x = [0 x_target x_target 0];
    y = [0 0 -1 -1];
    patch(x,y,'k','FaceAlpha',0.2)
    hold on
    plot(x_opt{1}(ind_t),x_opt{2}(ind_t),'linewidth',1.2,'color',0*ones(3,1));
    grid on
    hold on

    xlabel('$q_1$','interpreter','latex');
    ylabel('$q_2$','interpreter','latex');
    axis equal
    ylim([-0.4 max(x_opt{2})*1.15])
    xlim([0.0 x_target])
    saveas(gcf,'geometric_traj')
    %

    matlab_blue = [0 0.4470 0.7410];
    matlab_red = [0.8500 0.3250 0.0980];
    figure
    subplot(121)
    plot(x_opt{5},x_opt{3},'linewidth',1.2,'color',matlab_blue,'LineStyle','--');
    hold on
    plot(x_opt{5},x_opt{4},'linewidth',1.2,'color',matlab_red);

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
plot(x_opt_flat{1}(ind_t),x_opt_flat{2}(ind_t),'linewidth',1.2,'color',0*ones(3,1));
grid on
hold on
xlabel('$q_1$','interpreter','latex');
ylabel('$q_2$','interpreter','latex');
% axis equal
ylim([-0.4 max(x_opt_flat{2})*1.15])
xlim([0.0 x_target])
subplot(132)
plot(x_opt_flat{5},x_opt_flat{3},'linewidth',1.2,'color',matlab_blue);
hold on
plot(x_opt_flat{5},x_opt_flat{4},'linewidth',1.2,'color',matlab_red);
xlabel('$t$','interpreter','latex');
ylabel('$v(t)$','interpreter','latex');
grid on
legend({'$v_1(t)$','$v_2(t)$'},'interpreter','latex','NumColumns',2);
xlim([0 T]);
ylim([-6 6])
subplot(133)
stairs([x0(5),x_opt{5}(end, N_finite_elements:N_finite_elements:end)],[u1_opt;nan],'color',matlab_blue,'linewidth',1.2);
hold on
stairs([x0(5),x_opt{5}(end,N_finite_elements:N_finite_elements:end)],[u2_opt;nan],'color',matlab_red,'linewidth',1.2);
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

