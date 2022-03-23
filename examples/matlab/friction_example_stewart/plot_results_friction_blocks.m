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

for i = 1:n_u_DAE
    %     eval( ['mu_opt = alg_states(' num2str(i+2*n_theta) ':n_z+n_z*d:end);']);
    eval( ['u' num2str(i) '_opt = controls(' num2str(i) ':n_u:end);']);
end


if moving_finite_elements
    eval( ['h_opt = controls(n_u:n_u:end);']);
    tgrid = (cumsum([0;h_opt]));
    tgrid_z = cumsum(h_opt)';
end


practical_complementarity = [];
for i = 1:n_theta
    practical_complementarity = [practical_complementarity;eval( ['lambda' num2str(i) '_opt' '.*theta' num2str(i) '_opt ;'])];
end
fprintf('practical comp %4.2e \n',max(practical_complementarity))
%% controls
figure
stairs(tgrid_z,u1_opt)
hold on
plot(tgrid_z,F_input*cos(pi*tgrid_z))
xlabel('$t$','interpreter','latex');
ylabel('$u(t)$','interpreter','latex');
grid on

% optimal step size
if moving_finite_elements
    figure
    stairs(tgrid_z,h_opt)
    hold on
    stairs(tgrid_z,h_opt*0+h)
    xlabel('$t$','interpreter','latex');
    ylabel('$h_k(t)$','interpreter','latex');
    grid on
end
%% plots

%%
% algebraic

for ii = 1:n_simplex
    figure
    subplot(311)
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
    subplot(312)
    
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
    subplot(313)
     eval( ['plot(tgrid_z,mu' num2str(ii) '_opt);']);
    xlabel('$t$','interpreter','latex');
    ylabel(['$\mu_' num2str(ii) '(t)$' ],'interpreter','latex');
    grid on
end


%% Strong stationarity
figure
plot(theta1_opt+lambda1_opt)
hold on
plot(theta2_opt+lambda2_opt)
plot(theta3_opt+lambda3_opt)
grid on
yline(0,'k')
ylabel('$\theta+\lambda>0$','interpreter','latex');

a = min(theta1_opt+lambda1_opt);
b = min(theta2_opt+lambda2_opt);
c = min(theta3_opt+lambda3_opt);
biactivity = min([a;b;c])


%%
% figure
% stairs(h_eval(1,:)>h_eval(2,:))

figure
semilogy(complementarity_stats,'k')
xlabel('iter','interpreter','latex');
ylabel('complementarity','interpreter','latex');
grid on

%%
figure
subplot(611)
plot(tgrid_z,lambda1_opt);
hold on;
grid on
plot(tgrid_z,theta1_opt)
subplot(612)
stairs(tgrid_z,lambda1_opt.*theta1_opt);
grid on

subplot(613)
plot(tgrid_z,lambda2_opt);
hold on;
grid on
plot(tgrid_z,theta2_opt)
subplot(614)
stairs(tgrid_z,lambda2_opt.*theta2_opt);
grid on

subplot(615)
plot(tgrid_z,lambda3_opt);
hold on; 
plot(tgrid_z,theta3_opt)
grid on

subplot(616)
stairs(tgrid_z,lambda3_opt.*theta3_opt);
grid on
