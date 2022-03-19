clear all
clc
close all

import casadi.*
%% settings
% collocation settings
settings.d = 3;                            % Degree of interpolating polynomial
settings.collocation_scheme = 'radau';     % Collocation scheme: radau or legendre


% control stages


% MPCC settings
settings.solver_name = 'solver_fesd';           % name of solver function.
settings.mpcc_mode = 5;                    % 1 - extact, 2 - smooth  ,3 -relax , 4 - penalty, 5 - elastic mode
settings.s_elastic_max = 1e1;              % upper bound for elastic variables
settings.pointwise_or_integral = 1;        % 1 - pointwise relaxation of the comp., 0 integral relaxation of the cc.
settings.objective_scaling_direct = 1;                   % in penalty methods  1: J = J+(1/p)*J_comp (direct)  , 0 : J = p*J+J_comp (inverse)
% Penalty/Relaxation paraemetr
settings.sigma_0 = 1e0;                     % starting smouothing parameter
settings.sigma_N = 1e-10;                   % end smoothing parameter
settings.kappa = 0.1;                      % decrease rate
settings.N_homotopy = ceil(abs(log(settings.sigma_N/settings.sigma_0)/log(settings.kappa)))+1 ;% number of steps
settings.comp_tol = 1e-14;

% time freezing settings 
settings.time_freezing = 0;
settings.time_freezing_time_rescaling = 0;
settings.use_speed_of_time_variables =  0; % introduce s_tof for e
settings.local_speed_of_time_variable = 0;  
settings.stagewise_clock_constraint = 0;

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
%% Generate Model


model.active_control = 0; % OCP or simulation problem
model.use_hystereis_model = 1;
model.smooth_model = 0;
model.linear_auxiliary_dynamics = 1; 
model.time_optimal_problem = 0;
model.fuel_cost_on = 1;
model.fuel_cost_same = 1;

% model = car_hystheresis_model(model);
model = car_hystheresis_model_voronoi(model);
%% Reformulation of the PSS into a DCS
[model,settings] = model_reformulation_fesd(model,settings);

%% Call FESD Integrator - Simulation
model.T_sim = 3;
model.N_stages = 3;
model.N_finite_elements = 1;
model.h = 0.02;
settings.gamma_h = 1;
model.T = model.N_stages*model.h;
settings.use_previous_solution_as_initial_guess = 0;

[results,stats] = integrator_fesd(model,settings);

%% read data
unfold_struct(results,'base');
unfold_struct(stats,'base');
unfold_struct(model,'base');
tgrid = cumsum([0;h_vec]);
x1_opt =x_res(1,:);
x2_opt =x_res(2,:);
x3_opt =x_res(3,:);
x4_opt =x_res(4,:);
x5_opt =x_res(5,:);
%% complementarity stats
figure
% semilogy(tgrid(2:N_stages+1:end),complementarity_stats+1e-16,'k',LineWidth=1.5)
semilogy(complementarity_stats+1e-16,'k',LineWidth=1.5)
grid on
xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('Complementarity residual','Interpreter','latex')
%% in numerical time
figure
subplot(321)
plot(tgrid,x1_opt)
xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('$q(\tau)$','Interpreter','latex')
grid on
subplot(322)
plot(tgrid,x2_opt)
xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('$v(\tau)$','Interpreter','latex')
grid on
subplot(323)
plot(tgrid,x3_opt)
xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('$L(\tau)$','Interpreter','latex')
grid on
subplot(324)
plot(tgrid,x4_opt)
xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('$a(\tau)$','Interpreter','latex')
grid on
subplot(325)
plot(tgrid,x5_opt)
xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('$t(\tau)$','Interpreter','latex')
grid on


%% plots in phyisicaltime

figure
subplot(221)
plot(x5_opt,x1_opt)
xlabel('$t$ [phyisical time]','Interpreter','latex')
ylabel('$q(t)$','Interpreter','latex')
hold on
grid on
subplot(222)
plot(x5_opt,x2_opt)
grid on
yline(v1,'k--')
yline(v2,'k--')
xlabel('$t$ [phyisical time]','Interpreter','latex')
ylabel('$v(t)$','Interpreter','latex')
subplot(223)
plot(x5_opt,x3_opt)
grid on
xlabel('$t$ [phyisical time]','Interpreter','latex')
ylabel('$L(t)$','Interpreter','latex')

subplot(224)
plot(x5_opt,x4_opt)
grid on
xlabel('$t$ [phyisical time]','Interpreter','latex')
ylabel('$a(t)$','Interpreter','latex')
%% phase plot
figure
plot(x2_opt,x4_opt)
hold on
xline(v1,'k--')
xline(v2,'k--')
grid on
ylim([-0.1 1.1])
ylabel('$a$ [hystheresis state]','Interpreter','latex')
xlabel('$v$ [velocity]','Interpreter','latex')

%% plot with phase vectors
if model.use_hystereis_model
    plot_vector_fields_car_hysteresis
end