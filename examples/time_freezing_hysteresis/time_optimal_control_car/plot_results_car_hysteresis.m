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

function [ ] = plot_results_car_hysteresis(varargin)
close all
%[] = plot_results_throwing_ball(results,settings,model)
% TODO (anton.pozharskiy) this is still broken, come back and fix it before summer school.


%%
import casadi.*
if nargin == 0
    error('Results and model should be forwarded to this function.')
elseif nargin == 1
    results = varargin{1};
    unfold_struct(results,'caller');
elseif nargin == 2
    results = varargin{1};
    settings = varargin{2};
    unfold_struct(results,'caller');
    unfold_struct(settings,'caller');
elseif nargin == 3
    results = varargin{1};
    settings = varargin{2};
    model = varargin{3};
    unfold_struct(results,'caller');
    unfold_struct(model,'caller');
    unfold_struct(settings,'caller');
else
    results = varargin{1};
    settings = varargin{2};
    model = varargin{3};
    stats = varargin{4};
    unfold_struct(results,'caller');
    unfold_struct(model,'caller');
    unfold_struct(settings,'caller');
%     unfold_struct(stats,'caller');
end
dims = model.dims;
%%
N_finite_elements = settings.N_finite_elements(1);
h_k = h_k(1);
v1 = 10;
v2 = 15;
q_goal = 150;
v_goal = 0;
v_max = 25;
u_max = 5;

%%
%  Colors
extensive_plots = 0;
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
organe = [0.9290 0.6940 0.1250];
grey = [0.85 0.85 0.85];
%%
n = dims.n_x+dims.n_u;
nn = n-1;
%% Read solutions
if use_fesd
    h_opt = results.h;
    tgrid = results.t_grid;
end

%%
% figure
try
h_opt_stagewise = reshape(h_opt,N_finite_elements,N_stages);
if length(ind_tf)>1
    s_sot = results.w(ind_tf);
elseif length(ind_tf) == 0
    s_sot = ones(N_stages,1);
else
    s_sot = results.w(ind_tf)*ones(N_stages,1);
end
t_control_grid_pseudo = cumsum([0,sum(h_opt_stagewise)]);
t_control_grid_pseudo_streched = cumsum([0,sum(h_opt_stagewise).*s_sot']);
catch
    
end

%%
ind_t = find(diff([nan,results.x(5,:),nan])>1e-5);
time_physical_full = results.x(5,:);
time_physical = time_physical_full(ind_t);

%% plots in phyisical time for paper
figure
subplot(221)
plot(time_physical_full,results.x(2,:),'LineWidth',1.5)
hold on
plot(time_physical_full,results.x(2,:)*0+v1,'k--','LineWidth',1.0)
plot(time_physical_full,results.x(2,:)*0+v2,'k--','LineWidth',1.0)
plot(time_physical_full,results.x(2,:)*0+v_max,'r--','LineWidth',1.5)
xlabel('$t$ ','Interpreter','latex')
ylabel('$v(t)$ ','Interpreter','latex')
grid on

subplot(222)
stairs(time_physical_full(1:N_finite_elements:end),[results.u';nan],'LineWidth',1.5)
ylim([-u_max*1.1 u_max*1.1])
xlim([0 max(time_physical_full(1:N_finite_elements:end))])
xlabel('$t$ ','Interpreter','latex')
ylabel('$u(t)$ ','Interpreter','latex')
grid on

subplot(223)
plot(time_physical_full,results.x(4,:),'LineWidth',1.5)
xlabel('$t$ ','Interpreter','latex')
ylabel('$w(t)$ ','Interpreter','latex')
ylim([-0.1 1.1]);
grid on

ind_t = find(diff(time_physical_full)>0.01);
ind_t_complement = find(abs(diff(time_physical_full))<0.00000000000001);
x2_opt_phy = results.x(2,:);
x4_opt_phy = results.x(4,:);
x2_opt_phy(ind_t_complement) = nan;
x4_opt_phy(ind_t_complement) = nan;
subplot(224)
plot(results.x(2,:),results.x(4,:),'LineWidth',1.5)
hold on
% plot(x2_opt_phy,x4_opt_phy,linewidth=2)
% plot(x2_opt(ind_t),x4_opt(ind_t),linewidth=2)
grid on
ylim([-0.1 1.1])
ylabel('$w$','Interpreter','latex')
xlabel('$v$','Interpreter','latex')
tt = -3:0.5:3;
hold on
plot(tt*0+v1,tt,'k--')
plot(tt*0+v2,tt,'k--')
%     xline(v2)
saveas(gcf,'states_and_control')

%% phase plot

figure
    plot(results.x(2,:),results.x(4,:))
    grid on
    ylim([-0.1 1.1])
    ylabel('$w$ [hystheresis state]','Interpreter','latex')
    xlabel('$v$','Interpreter','latex')
    xline(v1)
    xline(v2)

%%  Homotopy complementarity stats
if nargin >= 4
figure
complementarity_stats = stats.complementarity_stats;
semilogy(complementarity_stats,'k','LineWidth',1.5)
xlabel('iter','interpreter','latex');
ylabel('Complementarity residual','interpreter','latex');
grid on
end

%% Plot homotopy results
if 0
    x_iter = sol.W(ind_x,:);
    % x_iter =x_iter(1:d+1:end,:);
    x_iter1= x_iter(1:dims.n_x:end,:);
    x_iter2= x_iter(2:dims.n_x:end,:);
    figure
    plot(x_iter1,x_iter2)

    legend_str = {};
    for ii = 1:size(x_iter2,2)
        legend_str  = [legend_str ; ['iter ' num2str(ii)]];
    end
    legend(legend_str);
end


end

