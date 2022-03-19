clear all
clc
close all

import casadi.*


%% settings
% collocation settings
settings = default_settings_fesd();

settings.d = 2;                            % Degree of interpolating polynomial
settings.collocation_scheme = 'radau';     % Collocation scheme: radau or legendre


% control stages

% MPCC settings
settings.solver_name = 'solver_fesd';           % name of solver function.
settings.mpcc_mode = 5;                    % 1 - extact, 2 - smooth  ,3 -relax , 4 - penalty, 5 - elastic mode
settings.s_elastic_max = 1e1;              % upper bound for elastic variables
settings.pointwise_or_integral = 1;        % 1 - pointwise relaxation of the comp., 0 integral relaxation of the cc.
settings.objective_scaling_direct = 0;                   % in penalty methods  1: J = J+(1/p)*J_comp (direct)  , 0 : J = p*J+J_comp (inverse)
% Penalty/Relaxation paraemetr
settings.sigma_0 = 1e1;                     % starting smouothing parameter
settings.sigma_N = 1e-10;                   % end smoothing parameter
settings.kappa = 0.1;                      % decrease rate
settings.N_homotopy = ceil(abs(log(settings.sigma_N/settings.sigma_0)/log(settings.kappa)))+1 ;% number of steps
settings.comp_tol = 1e-14;

% time freezing settings 
settings.time_freezing = 0;
settings.time_optimal_problem = 0;
settings.time_rescaling = 0;
settings.use_speed_of_time_variables =  0; % introduce s_tof for e
settings.local_speed_of_time_variable = 0;  
settings.stagewise_clock_constraint = 0;
settings.initial_theta = 0.5;

%^ IPOPT Settings
opts_ipopt.verbose = false;
opts_ipopt.ipopt.max_iter = 250;
% opts_ipopt.ipopt.max_iter = 100;
opts_ipopt.ipopt.print_level = 0;
% opts_ipopt.ipopt.print_level = 5;
opts_ipopt.ipopt.mu_strategy = 'adaptive';
opts_ipopt.ipopt.mu_oracle = 'quality-function';
settings.opts_ipopt = opts_ipopt;
% opts_ipopt.ipopt.tol = 1e-5;
% 1 - scale inverse, 0 scale direct


% finite elements with switch detection
settings.use_fesd = 1;       % turn on moving finite elements algortihm
settings.gamma_h = 1;                    % how much can h addapt
settings.regularize_h = 0;                 % add quadratic cost term rho_h = (h-hk)^2
settings.rho_h = 1;                        % regularization penalty
settings.piecewise_equidistant_grid = 1;
settings.fesd_complementartiy_mode = 1;
settings.delta_h_regularization = 0;
settings.piecewise_equidistant_grid_sigma = 0.1;
settings.equidistant_control_grid = 0;

%% Time settings


%% Generate Model
% model = temp_control_model();
model = temp_control_model_voronoi();

%% Reformulation of the PSS into a DCS
[model,settings] = model_reformulation_fesd(model,settings);

%% Call FESD Integrator - Simulation
model.T_sim = 4;
model.N_stages = 2;
model.N_finite_elements = 1;
model.h = 0.01;
settings.gamma_h = 1;
model.T = model.N_stages*model.h;
settings.use_previous_solution_as_initial_guess = 1;
[results,stats] = integrator_fesd(model,settings);

%% read data
unfold_struct(results,'base');
unfold_struct(stats,'base');
unfold_struct(model,'base');
tgrid = cumsum([0;h_vec]);
x1_opt =x_res(1,:);
x2_opt =x_res(2,:);
x3_opt =x_res(3,:);
%% complementarity stats
figure
% semilogy(tgrid(2:N_stages+1:end),complementarity_stats+1e-16,'k',LineWidth=1.5)
semilogy(complementarity_stats+1e-16,'k',LineWidth=1.5)
grid on
xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('Complementarity residual','Interpreter','latex')
%%

%% in numerical time
t_phy = [x3_opt,nan];
dt_phy = diff(t_phy);
ddt_phy = diff([dt_phy,nan]);
ind_t = find(abs(dt_phy)>0.01);
ind_t_complement = find(abs(dt_phy)<=0.001);
ind_dt = find(abs(ddt_phy)<0.000005);

time_frozen = tgrid*0;
time_frozen(ind_dt)=5;

% junction point
figure
subplot(311)
plot(tgrid,x1_opt)
% xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('$x(\tau)$','Interpreter','latex')
grid on
subplot(312)
plot(tgrid,x2_opt)
% xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('$w(\tau)$','Interpreter','latex')
grid on
subplot(313)
plot(tgrid,x3_opt)
hold on
area(tgrid,time_frozen,FaceAlpha=0.1)
ylim([0 2])
xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('$t(\tau)$ [phyisical time]','Interpreter','latex')
grid on
% phase plot
figure
x1_phy = x1_opt;
x1_phy(ind_t_complement) = nan;
x2_phy = x2_opt;
x2_phy(ind_t_complement) = nan;
plot(x1_opt,x2_opt)
hold on
plot(x1_phy,x2_phy,'LineWidth',1.5)
xline(y1,'k--')
xline(y2,'k--')
grid on
ylim([-0.1 1.1])
ylabel('$w$','Interpreter','latex')
xlabel('$x$ ','Interpreter','latex')
% plots in phyisicaltime
figure
subplot(211)
plot(x3_opt,x1_opt)
hold on
xlabel('$t$ [phyisical time]','Interpreter','latex')
ylabel('$\theta(t)$ [temperature]','Interpreter','latex')
hold on
grid on
yline(y1,'k--')
yline(y2,'k--')
subplot(212)
plot(x3_opt,x2_opt)
grid on
%%
figure
for ii = 1:4
subplot(4,1,ii)    
stairs(tgrid,[nan,theta_res(ii,:)]);
grid on
ylabel(['$\theta_' num2str(ii) '(t)$'],'Interpreter','latex')
ylim([-0.2 1.2]);
end
xlabel('$\tau$','Interpreter','latex')

%% Read and plot Result
% plot_results_for_paper
% plot_results_temp_control

ind_res  = [];
g_res = [];
for ii  = 1:length(x_res)
ind_res = [ind_res ;find(full(h_fun(x_res(:,ii))) <= full(min(h_fun(x_res(:,ii)))))];
g_res  = [g_res,full(h_fun(x_res(:,ii)))];
end
figure
plot(tgrid,g_res)