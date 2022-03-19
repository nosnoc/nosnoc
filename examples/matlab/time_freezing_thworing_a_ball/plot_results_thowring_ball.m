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
%% read solutions

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
    
%     eval( ['h_opt = controls(n_u:n_u:end);']);
    h_opt = w_opt(ind_h);
    tgrid = (cumsum([0;h_opt]));
    tgrid_z = cumsum(h_opt)';
end

%%
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

%%

%%
% figure
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
figure
stairs(tgrid(1:N_finite_elements:end),[s_sot;nan],'linewidth',1.2)
xlabel('$\tau$','interpreter','latex');
ylabel('$s(\tau)$','interpreter','latex');
grid on
ylim([0.5 3.2])

%% 
if extensive_plots
figure
ax(1) = subplot(211);
plot(tgrid,x5_opt*2,'linewidth',2)
hold on

stairs(t_control_grid_pseudo,[s_sot;nan],'linewidth',2)
yline(1,'k')
xlabel('$\tau$','interpreter','latex');
for ii = 2:N_stages+1
    xline(t_control_grid_pseudo(ii),'k--');
end
legend({'$2 t(\tau)$','$s_{\mathrm{sot}}(\tau)$','$s_{\mathrm{sot}} = 1$'},'interpreter','latex');
ax(2) = subplot(212);
stairs(t_control_grid_pseudo,[u1_opt;nan],'linewidth',2)
hold on
stairs(t_control_grid_pseudo,[u2_opt;nan],'linewidth',2)
xlabel('$\tau$','interpreter','latex');
ylabel('$u(\tau)$','interpreter','latex')
for ii = 2:N_stages+1
    xline(t_control_grid_pseudo(ii),'k--');
end
linkaxes(ax,'x')
% set(ax,'xlim',[0 1])

end

%% with steched time
if 0 
figure

temp = h_opt_stagewise.*s_sot';
temp = temp(:);
tgrid_streched = cumsum([0; temp]);


ax(1) = subplot(211);
plot(tgrid_streched,x5_opt*2,'linewidth',2)
hold on
stairs(t_control_grid_pseudo_streched,[s_sot;nan],'linewidth',2)
yline(1,'k')
xlabel('$\tau$','interpreter','latex');
for ii = 2:N_stages+1
    xline(t_control_grid_pseudo_streched(ii),'k--');
end
legend({'$2 t(\tau)$','$s_{\mathrm{sot}}(\tau)$','$s_{\mathrm{sot}} = 1$'},'interpreter','latex');
ax(2) = subplot(212);
stairs(t_control_grid_pseudo_streched,[u1_opt;nan],'linewidth',2)
hold on
stairs(t_control_grid_pseudo_streched,[u2_opt;nan],'linewidth',2)
xlabel('$\tau$','interpreter','latex');
ylabel('$u(\tau)$','interpreter','latex')
for ii = 2:N_stages+1
    xline(t_control_grid_pseudo_streched(ii),'k--');
end
linkaxes(ax,'x')
% set(ax,'xlim',[0 1])

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
x = [0 6 6 0];
y = [0 0 -1 -1];
patch(x,y,'red','FaceAlpha',0.2)
hold on
plot(x1_opt(ind_t),x2_opt(ind_t),'linewidth',1.2,'color',0*ones(3,1));
grid on
hold on
% ttt = 0:0.2:6;
% ha = area(ttt,0*ttt-1);
% ha.FaceAlpha = 0.5;

xlabel('$q_x$','interpreter','latex');
ylabel('$q_y$','interpreter','latex');
axis equal
ylim([-0.4 1.5])
xlim([0.0 6])
%%
figure
plot(x5_opt,x3_opt,'linewidth',1.2,'color',0*ones(3,1),'LineStyle','--');
hold on
plot(x5_opt,x4_opt,'linewidth',1.2,'color',0*ones(3,1));
xlabel('$t$','interpreter','latex');
ylabel('$v(t)$','interpreter','latex');
grid on
legend({'$v_x(t)$','$v_y(t)$'},'interpreter','latex');
%%
figure
stairs(x5_opt(1:N_finite_elements:end),[u1_opt;nan],'color',0*ones(3,1),'linewidth',1.2);
hold on
stairs(x5_opt(1:N_finite_elements:end),[u2_opt;nan],'color',0*ones(3,1),'LineStyle','--','linewidth',1.2);
grid on
xlabel('$t$','interpreter','latex');
ylabel('$u(t)$','interpreter','latex');
grid on
legend({'$u_x(t)$','$u_y(t)$'},'interpreter','latex');
%



%% Algebraic variavbles
%%
% make_animation_and_video
%%

for lll= 1:1; 

if 0
    %%
    
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



%% all in one plot
if 0 
streched_time = 0;
if streched_time
    t_control_grid_pseudo = t_control_grid_pseudo_streched;
end
figure
ax(1) = subplot(611);
stairs(t_control_grid_pseudo,[s_sot;nan])
hold on
yline(1,'r')
hold on
xlabel('$\tau$','interpreter','latex');
ylabel('$s_{\mathrm{sot}}(\tau)$','interpreter','latex');
for ii = 2:N_stages+1
    xline(t_control_grid_pseudo(ii),'k--');
end
ax(2) = subplot(612);
% step sizes
stairs(tgrid,[h_opt;nan]);
xlabel('$\tau$','interpreter','latex');
ylabel('$h_{k,i}$','interpreter','latex');
for ii = 2:N_stages+1
    xline(t_control_grid_pseudo(ii),'k--');
end
ax(3) = subplot(613);
% states
plot(tgrid,x1_opt);
hold on
plot(tgrid,x2_opt);
xlabel('$\tau$','interpreter','latex');
ylabel('$q(\tau)$','interpreter','latex');
for ii = 2:N_stages+1
    xline(t_control_grid_pseudo(ii),'k--');
end
ax(4) = subplot(614);
% velocity
plot(tgrid,x3_opt);
hold on
plot(tgrid,x4_opt);
xlabel('$\tau$','interpreter','latex');
ylabel('$v(\tau)$','interpreter','latex')
for ii = 2:N_stages+1
    xline(t_control_grid_pseudo(ii),'k--');
end
ax(5) = subplot(615);
% controls
stairs(t_control_grid_pseudo,[u1_opt;nan])
hold on
stairs(t_control_grid_pseudo,[u2_opt;nan])
xlabel('$\tau$','interpreter','latex');
ylabel('$u(\tau)$','interpreter','latex')
for ii = 2:N_stages+1
    xline(t_control_grid_pseudo(ii),'k--');
end
ax(6) = subplot(616);
% clock state
 plot(tgrid,x5_opt);
xlabel('$\tau$','interpreter','latex');
ylabel('$t(\tau)$','interpreter','latex');   
for ii = 2:N_stages+1
    xline(t_control_grid_pseudo(ii),'k--');
end

linkaxes(ax,'x')
set(ax,'xlim',[0 2])


end
end