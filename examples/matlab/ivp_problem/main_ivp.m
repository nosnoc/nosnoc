clear all
clc
close all
import casadi.*
%%
objective_fun_experiment = 0;
local_minima_experiment = 1;
N_samples = 100;
%% settings
settings = default_settings_fesd();
settings.n_s = 2;     
settings.s_elastic_0 = 1e1;
settings.mpcc_mode = 5;           
settings.sigma_0 = 1e-4;
settings.sigma_N = 1e-5;
settings.N_homotopy = 1;
% settings.irk_scheme = 'legendre';
settings.cross_comp_mode = 8;
settings.pss_mode = 'Step';
settings.use_fesd = 0;
settings.equidistant_control_grid = 0;
%% Generate Model
model = ivp_problem();
%% Solve NLP
%  [results,stats,model,settings,solver_initalization] = nosnoc_solver(model,settings);
% complementarity_stats = stats.complementarity_stats;
% fprintf('Total CPU Time = %2.4f : \n',sum(stats.cpu_time));
%% objective expereimnt
if objective_fun_experiment
    objective = [];
    complementarity_stats_obj = [];
    x0_vec = linspace(-2,-0.8,N_samples);
    error_state = [];
    error_objective = [];
    [solver,solver_initalization, model,settings] = create_nlp_fesd(model,settings);
    unfold_struct(model,'base');
    unfold_struct(settings,'base');
    unfold_struct(solver_initalization,'base');

    for jj = 1:N_samples
        x0 = x0_vec(jj);
        solver_initalization.lbw(1) =  x0;
        solver_initalization.ubw(1) =  x0;
        % solve NLP
        [sol,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);
        w_opt = full(sol.x);
        x1_opt = w_opt(ind_x);
        solver_initalization.w0 = w_opt;
        f_opt = full(J_fun(w_opt));
        objective = [objective;f_opt];
        complementarity_iter = full(comp_res(w_opt));
        complementarity_stats_obj = [complementarity_stats_obj;complementarity_iter ];
        ts = -x0/3;
        L_analytic =(8/3*ts^3 +8/3*x0*ts^2 + 8/9*x0^2*ts +1/3*T^3+1/3*x0*T^2+x0^2*T/9) + (T+(x0-5)/3)^2;
        L_analytic =(8/3*ts^3 +8/3*x0*ts^2 + 8/9*x0^2*ts +1/3*T^3+1/3*x0*T^2+x0^2*T/9) + (T+(x0-5)/3)^2;
        error_state = [error_state;abs(x1_opt(end)-(T-ts)) ];
        error_objective = [error_objective;abs(f_opt-L_analytic)];
    end

    figure
    subplot(211)
    semilogy(x0_vec,error_state)
    grid on
    ylabel('Error state','interpreter','latex');
    subplot(212)
    semilogy(x0_vec,error_objective)
    grid on
    xlabel('$x_0$','interpreter','latex');
    ylabel('Error objective','interpreter','latex');
    %
    figure
    plot(x0_vec,objective,'LineWidth',1.0);
    ts = -x0_vec/3;
    L_analytic =(8/3*ts.^3 +8/3*x0_vec.*ts.^2 + 8/9*x0_vec.^2.*ts +1/3*T^3+1/3*x0_vec.*T^2+x0_vec.^2*T/9) + (T+(x0_vec-5)/3).^2;
    hold on
    plot(x0_vec,L_analytic,':','LineWidth',2.5);
    xlabel('$x_0$','interpreter','latex');
    ylabel('Objective','interpreter','latex');
    if settings.use_fesd
    legend({'$V_{\mathrm{FESD}}(x_0)$','$V^{*}(x_0)$'},'interpreter','latex');
    else
      legend({'$V_{\mathrm{Num}}(x_0)$','$V^{*}(x_0)$'},'interpreter','latex');
    end
    grid on
end

%% Convergence test:
if local_minima_experiment 
    objective = [];
    complementarity_stats_obj = [];
    x0_vec = linspace(-2,-0.8,N_samples);

    [solver,solver_initalization, model,settings] = create_nlp_fesd(model,settings);
    unfold_struct(model,'base');
    unfold_struct(settings,'base');
    unfold_struct(solver_initalization,'base');
    
    x0_star = [];
    solver_initalization.lbw(1) =  -inf;
    solver_initalization.ubw(1) =  inf;

    for jj = 1:N_samples
        x0 = x0_vec(jj);
        solver_initalization.w0(ind_x) = x0;
        solver_initalization.lbw(1) =  x0;
        solver_initalization.ubw(1) =  x0;
        %initalizer
        [sol,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);

        % solve NLP
        solver_initalization.w0 = full(sol.x);
        solver_initalization.lbw(1) =  -inf;
        solver_initalization.ubw(1) =  inf;
        [sol,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);
        w_opt = full(sol.x);
        x1_opt = w_opt(ind_x);

        x0_star = [x0_star;x1_opt(1)];
        f_opt = full(J_fun(w_opt));
        objective = [objective;f_opt];
        complementarity_iter = full(comp_res(w_opt));
        complementarity_stats_obj = [complementarity_stats_obj;complementarity_iter];
        ts = -x0/3;
        L_analytic =(8/3*ts^3 +8/3*x0*ts^2 + 8/9*x0^2*ts +1/3*T^3+1/3*x0*T^2+x0^2*T/9) + (T+(x0-5)/3)^2;
%         error_state = [error_state;abs(x1_opt(end)-(T-ts)) ];
%         error_objective = [error_objective;abs(f_opt-L_analytic)];
    end

%     x0_analytic = -1.427570148502220;
    x0_analytic = -1.427572232082767;
%      = double(x_star(1))
    figure
    plot(x0_vec,x0_analytic*ones(1,N_samples),'LineWidth',1.0);
    hold on
    plot(x0_vec,x0_star,'LineWidth',1.0);
    axis equal
    xlabel('$x_0$','interpreter','latex');
    ylabel('$x_0^*$','interpreter','latex');
    if settings.use_fesd
        legend({'Analytic Solution','FESD'},'interpreter','latex');
    else
        legend({'Analytic Solution','Numerical Solution'},'interpreter','latex');
    end
    grid on
end

%% Read and plot Result
% plot_results_ivp

