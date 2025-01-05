%%
clear all; clc; close all;
import casadi.*
%% basic settings
plot_integrator_output = 1;
plot_continious_time_sol = 1;
smooth_model = 0; % if 1, use model without switch for a sanity check
%% discretization settings
T_sim = pi/2;
N_sim  = 29;
N_FE = 2;

%% Init
problem_options = nosnoc.Options();
integrator_opts = nosnoc.integrator.Options();
solver_options = integrator_opts.fesd_solver_opts; % the fesd integrator uses an mpec solver, call and modify its options
model = nosnoc.model.Pss();


% select integrator
% integrator_options.integrator_plugin = "FESD";
% integrator_options.integrator_plugin = "SMOOTHED_PSS";
integrator_opts.integrator_plugin = "STEWART";

% integraotr options
problem_options.dcs_mode = 'Stewart'; % 'Step;
% problem_options.use_fesd = 1;       % switch detection method on/off
% problem_options.rk_scheme = RKSchemes.RADAU_IIA; %'Gauss-Legendre';
% problem_options.n_s = 4;
problem_options.N_finite_elements = N_FE; % number of finite elements
problem_options.T_sim = T_sim; % total simulation times
problem_options.N_sim = N_sim; % number of simulations step during T_sim

% MPEC solver options
solver_options.print_level = 2;
solver_options.homotopy_steering_strategy = HomotopySteeringStrategy.ELL_INF;
solver_options.complementarity_tol = 1e-9; % Penalty/Relaxation parameter

%% Set up model;
x0 = [exp(-1);0]; % inital value
x_star = [exp(T_sim-1)*cos(2*pi*(T_sim-1));-exp((T_sim-1))*sin(2*pi*(T_sim-1))]; % Analytic solution, if T > 1;
% model parameters
omega = -2*pi;
R_osc  = 1;

A1 = [1 omega;...
    -omega 1];
A2 = [1 -omega;...
    omega 1];
if smooth_model
    A2 = A1;
end

% Variable defintion
x = SX.sym('x',2);
c = x'*x-1; % the switching surface is the unit circle
f_1 = A1*x;
f_2 = A2*x;
F = [f_1 f_2];

% Populate model
model.F = F;
model.x = x;
model.c = c;
model.S = [-1;1];
model.x0 = x0;

%% Create Stewart problem manually

% fill in model
model.verify_and_backfill(problem_options);

% compute current epsilon active set;
% min(g(xk)) and use epsilion from options;
% I_k % current active set
I_active = [1];
I_all = 1:4;
I_inactive = [];
f_x_stewart = 0;
ii = 1;
% for ii = 1:model.dims.n_sys
% z_full = zeros(model.dims.n_x);
F_ii = model.F{ii};
g_indicator_ii = model.g_indicator{ii};
nabla_g = g_indicator_ii.jacobian(model.x);
M_full = nabla_g*model.F{ii};
alpha = SX.sym('alpha');
M_alpha = M_full+alpha;
M_alpha_fun = Function('M_alpha_fun', {model.x, model.u, alpha}, {M_alpha});
z_hat = inv(M_alpha(I_k,I_k))*ones(length(I_k),1);
z_hat = z_hat./sum(z_hat);

% z_full -
f_x_stewart = f_x_stewart+F_ii(:,I_k)*z_hat;
rhs_fun = Function('rhs_fun', {model.x, model.u, alpha}, {f_x_stewart});
ode_func = @(t, x , u ,alpha) full(rhs_fun(x, u, alpha));

g_active = g_indicator_ii(I_k);
g_inactive = g_indicator_ii;
g_inactive(I_k) = [];

psi_switching_fun = min([min(g_active)-min(g_inactive);z_hat]);
switching_fun = Function('switching_fun', {model.x, model.u, alpha}, {psi_switching_fun});
switching_fun = @(t, x , u , alpha) full(switching_fun(x, u, alpha));

% esimate alpha
% create M_alpha matrix only for I
% write formula for z reduced; and zeros for others
% create ode func;
% rhs_fun = Function('rhs_fun', {model.x, model.u, sigma}, {f_x_smoothed});
%                 obj.ode_func = @(t, x , u ,sigma) full(rhs_fun(x, u, sigma));
% LCP

% end


obj = struct;
x0 = model.x0;
obj.x_res = x0;
obj.x_res_full = x0;
obj.t_grid = 0;
obj.t_grid_full = 0;
obj.x_curr = x0;
obj.ode_func = ode_func;
% obj.set_x0(x0);
t_current = 0;
T = T_sim/N_sim;

for ii=1:N_sim
    u_i = [];
    switch(integrator_opts.matlab_ode_solver)
        case 'ode45'
            [t_sim, x_sim] = ode45(@(t, x)  obj.ode_func(t,x,u_i,integrator_opts.alpha_stewart), [t_current, t_current+opts.T], obj.x_curr, integrator_opts.matlab_ode_opts);
        case 'ode23'
            [t_sim, x_sim] = ode23(@(t, x)  obj.ode_func(t,x,u_i,integrator_opts.alpha_stewart), [t_current, t_current+opts.T], obj.x_curr, integrator_opts.matlab_ode_opts);
        case 'ode113'
            [t_sim, x_sim] = ode113(@(t, x)  obj.ode_func(t,x,u_i,integrator_opts.alpha_stewart), [t_current, t_current+opts.T], obj.x_curr, integrator_opts.matlab_ode_opts);
        case 'ode78'
            [t_sim, x_sim] = ode78(@(t, x)  obj.ode_func(t,x,u_i,integrator_opts.alpha_stewart), [t_current, t_current+opts.T], obj.x_curr, integrator_opts.matlab_ode_opts);
        case 'ode89'
            [t_sim, x_sim] = ode89(@(t, x)  obj.ode_func(t,x,u_i,integrator_opts.alpha_stewart), [t_current, t_current+opts.T], obj.x_curr, integrator_opts.matlab_ode_opts);
        case 'ode15s'
            [t_sim, x_sim] = ode15s(@(t, x)  obj.ode_func(t,x,u_i,integrator_opts.alpha_stewart), [t_current, t_current+opts.T], obj.x_curr, integrator_opts.matlab_ode_opts);
        case 'ode23s'
            [t_sim, x_sim] = ode23s(@(t, x)  obj.ode_func(t,x,u_i,integrator_opts.alpha_stewart), [t_current, t_current+T], obj.x_curr, integrator_opts.matlab_ode_opts);
        case 'ode23t'
            [t_sim, x_sim] = ode23t(@(t, x)  obj.ode_func(t,x,u_i,integrator_opts.alpha_stewart), [t_current, t_current+opts.T], obj.x_curr, integrator_opts.matlab_ode_opts);
        case 'ode23tb'
            [t_sim, x_sim] = ode23tb(@(t, x)  obj.ode_func(t,x,u_i,integrator_opts.alpha_stewart), [t_current, t_current+opts.T], obj.x_curr, integrator_opts.matlab_ode_opts);
        case {'cvodesnonstiff','cvodesstiff','idas'} % Handle with new OO ode method.\
            % Note: these are only available >=2024a
            % TODO(@anton) maybe instead of looping use Refine=N.
            if ~exist('F', 'var')
                F = ode;
                F.ODEFcn = @(t, x)  obj.ode_func(t,x,u_i,integrator_opts.alpha_stewart);
                F.Solver = obj.integrator_opts.matlab_ode_solver;
                F.AbsoluteTolerance = odeget(integrator_opts.matlab_ode_opts, 'AbsTol', 1e-6);
                F.RelativeTolerance = odeget(integrator_opts.matlab_ode_opts, 'RelTol', 1e-3);
            end
            F.InitialTime = t_current;
            F.InitialValue = obj.x_curr';

            sol = F.solve(t_current, t_current+opts.T);
            t_sim = sol.Time';
            x_sim = sol.Solution';
    end
    t_current = t_sim(end);
    obj.t_grid = [obj.t_grid, t_sim(end)];
    obj.t_grid_full = [obj.t_grid_full t_sim(2:end)'];
    obj.x_res_full = [obj.x_res_full, x_sim(2:end,:)'];
    obj.x_res = [obj.x_res, x_sim(end,:)'];
    % obj.set_x0(x_sim(end,:)');
    obj.x0 = (x_sim(end,:)');
    obj.x_curr = obj.x0;
end
x_res = obj.x_res;
x_res_full = obj.x_res_full;
t_grid = obj.t_grid;
t_grid_full = obj.t_grid_full;

%% Call integrator
% integrator = nosnoc.Integrator(model, problem_options, integrator_options);
% [t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();

%% numerical error
x_fesd = x_res(:,end);
error_x = norm(x_fesd-x_star,"inf");
fprintf(['Numerical error with h = %2.3f and ' char(problem_options.rk_scheme) ' with n_s = %d stages is: %5.2e: \n'],problem_options.h_sim,problem_options.n_s,error_x);
%% plot_solution_trajectory
t_star = R_osc; % eact switching time

h_opt_full = diff(obj.t_grid);

x1_opt = x_res(1,:);
x2_opt = x_res(2,:);
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

