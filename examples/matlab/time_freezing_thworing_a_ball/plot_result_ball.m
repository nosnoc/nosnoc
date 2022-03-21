function [varargout] = plot_result_ball(model,settings,results,stats)

unfold_struct(model,'caller');
unfold_struct(settings,'caller')
d = 3;
% unfold_struct(resultsver_initalization,'caller');

obj = full(results.f);
w_opt = full(results.x);

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
%% read resultsutions

diff_states = w_opt(ind_x);
controls = w_opt(ind_u);
alg_states = w_opt(ind_z);

% differential states
for i = 1:n_x
    eval( ['x' num2str(i) '_opt = diff_states(' num2str(i) ':n_x+n_x*d:end);']);
end
%
% for i = 1:n_x
%     eval( ['x' num2str(i) '_opt = diff_states(' num2str(i) ':n_x:end);']);
% end

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
% lambdas
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

time_physical = x5_opt(ind_t);
% Geomtric plot

figure;
x = [0 4 4 0];
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
xlim([0.0 4])
%%
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



%% Algebraic variavbles

if 0

    % algebraic

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



    %%
    complementarity_stats = stats.complementarity_stats;
    figure
    semilogy(complementarity_stats,'k')
    xlabel('iter','interpreter','latex');
    ylabel('complementarity','interpreter','latex');
    grid on
end

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
    % xlim([0 0.6])
    % ylim([0 0.6])
    axis equal
    xlim([0 T])
    ylim([0 T])
    xlabel('$\tau$ [Numerical Time]','interpreter','latex');
    ylabel('$t(\tau)$ [Phyisical Time]','interpreter','latex');
end
end

