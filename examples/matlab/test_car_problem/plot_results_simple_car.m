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
function [ ] = plot_results_simple_car(varargin)
close all
%[] = plot_results_throwing_ball(results,settings,model)
% crate and solve and OCP for an example time freezing problem
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

%% TODO: improve ploting codes for variable finte elements.
N_finite_elements = N_finite_elements(1); 
h_k = h_k(1);

%%
%  Colors
extensive_plots = 0;
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
organe = [0.9290 0.6940 0.1250];
grey = [0.85 0.85 0.85];
%%
n = n_x+n_u;
nn = n-1;
tgrid = linspace(0, T, N_stages+1);
tgrid_z = linspace(0, T, N_stages);
%% Read solutions

diff_states = w_opt(ind_x);
controls = w_opt(ind_u);
alg_states = w_opt(ind_z);

% differential states
for i = 1:n_x
    eval( ['x' num2str(i) '_opt = diff_states(' num2str(i) ':n_x+n_x*n_s:end);']);
end

% convex multiplers
for i = 1:n_theta
    eval( ['theta' num2str(i) '_opt = alg_states(' num2str(i) ':n_z+n_z*(n_s-1):end);']);
end
% lambdas
for i = 1:n_theta
    eval( ['lambda' num2str(i) '_opt = alg_states(' num2str(i+n_theta) ':n_z+n_z*(n_s-1):end);']);
end
% mu
for i = 1:n_simplex
    eval( ['mu' num2str(i) '_opt = alg_states(' num2str(i+2*n_theta) ':n_z+n_z*(n_s-1):end);']);
end

for i = 1:n_u
    eval( ['u' num2str(i) '_opt = controls(' num2str(i) ':n_u:end);']);
end


if use_fesd

    %     eval( ['h_opt = controls(n_u:n_u:end);']);
    h_opt = w_opt(ind_h);
    tgrid = (cumsum([0;h_opt]));
    tgrid_z = cumsum(h_opt)';
end


%%
% figure
if 1

    h_opt_stagewise = reshape(h_opt,N_finite_elements,N_stages);
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


%% plots in phyisicaltime
figure
subplot(221)
plot(x5_opt,x1_opt)
xlabel('$t$ ','Interpreter','latex')
ylabel('$q(t)$ ','Interpreter','latex')
grid on
subplot(222)
plot(x5_opt,x2_opt)
hold on
plot(x5_opt,x2_opt*0+v_trash_hold,'k-')
xlabel('$t$ ','Interpreter','latex')
ylabel('$v(t)$ ','Interpreter','latex')

grid on

subplot(223)
plot(x5_opt,x4_opt)
xlabel('$t$ ','Interpreter','latex')
ylabel('$w(t)$ ','Interpreter','latex')
ylim([-0.1 1.1]);
grid on
subplot(224)
stairs(x5_opt(1:N_finite_elements:end),[u1_opt;nan])
xlabel('$t$ ','Interpreter','latex')
ylabel('$u(t)$ ','Interpreter','latex')
grid on



%%
if mpcc_mode == 4
    ind_t = find([1;theta1_opt]>1e-2);
else
    ind_t = find(diff([nan;x5_opt;nan])>1e-5);
end

time_physical = x5_opt(ind_t);




%% Plot homotopy results.
if 0
    x_iter = sol.W(ind_x,:);
    % x_iter =x_iter(1:d+1:end,:);
    x_iter1= x_iter(1:n_x:end,:);
    x_iter2= x_iter(2:n_x:end,:);
    figure
    plot(x_iter1,x_iter2)

    legend_str = {};
    for ii = 1:size(x_iter2,2)
        legend_str  = [legend_str ; ['iter ' num2str(ii)]];
    end
    legend(legend_str);
end

%% Algebraic variavbles
if 1
    for ii = 1:n_simplex
        figure
        subplot(311);
        for i = m_ind_vec(ii):m_ind_vec(ii)+m_vec(ii)-1
            eval( ['plot(tgrid_z,theta' num2str(i) '_opt);']);
            if  i == m_ind_vec(ii)
                hold on
            end
        end
        xlabel('$t$','interpreter','latex');
        ylabel('$\theta(t)$','interpreter','latex');
        hold on
        grid on
        legend_str= {};
        for i = m_ind_vec(ii):m_ind_vec(ii)+m_vec(ii)-1
            legend_str = [legend_str, ['$\theta_' num2str(i) '(t)$']];
        end
        legend(legend_str ,'interpreter','latex');
        subplot(312);
        for i = m_ind_vec(ii):m_ind_vec(ii)+m_vec(ii)-1
            eval( ['plot(tgrid_z,lambda' num2str(i) '_opt);']);
            if  i == m_ind_vec(ii)
                hold on
            end
        end
        xlabel('$t$','interpreter','latex');
        ylabel('$\lambda(t)$','interpreter','latex');
        hold on
        grid on
        legend_str= {};
        for i = m_ind_vec(ii):m_ind_vec(ii)+m_vec(ii)-1
            legend_str = [legend_str, ['$\lambda' num2str(i) '(t)$']];
        end
        legend(legend_str ,'interpreter','latex');
        subplot(313);
        eval( ['plot(tgrid_z,mu' num2str(ii) '_opt);']);
        xlabel('$t$','interpreter','latex');
        ylabel(['$\mu_' num2str(ii) '(t)$' ],'interpreter','latex');
        grid on
    end

end


end

