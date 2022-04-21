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
clear all
clc
close all
import casadi.*
%%
model = [];
[settings] = default_settings_fesd();  %% Optionally call this function to have an overview of all options.
settings.print_level = 3;
settings.mpcc_mode = 3;
settings.N_homotopy = 15;
settings.comp_tol = 1e-12;
settings.irk_representation = 'differential';
N_stages = 6;
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 2;
scenario.N_finite_elements = 2;

settings.use_fesd = 1;

scenario.save_results = 1;
scenario.estimate_terminal_error = 1;


scenario.illustrate_regions  = 1;
scenario.terminal_constraint = 1;
scenario.linear_control = 1;
scenario.rho_v = 1;
scenario.rho_u = 1;
scenario.plot_results_sliding_mode = 0;


%% Basic
mpcc_mode_fesd = [5 5 5 3 3];
mpcc_mode_std = [3 3 3 3 3];
irk_schemes = {'Radau-IIA','Lobatto-IIIA','Lobatto-IIIC','Gauss-Legendre','Explicit-RK'};
n_s_vec = {[1:4],[2:4],[2:4],[1:4],[1 4]};
orders_vec = {n_s_vec{1}*2-1,n_s_vec{2}*2-2,n_s_vec{2}*2-2,n_s_vec{3}*2,n_s_vec{4}};
N_finite_elements_vec = 2:2:12;
N_experiments = 1:5;

%%
N_experiments = [1 2 4];
N_finite_elements_vec = 2:3;
n_s_vec = {[1:4],[2:4],[2:3]};
%%
mpcc_mode_fesd = [5 5 3 3];
mpcc_mode_std = [3 3 3 3];
irk_schemes = {'Radau-IIA','Lobatto-IIIC','Gauss-Legendre','Explicit-RK'};
n_s_vec = {[1:4],[2:4],[1:4],[1 4]};
orders_vec = {n_s_vec{1}*2-1,n_s_vec{2}*2-2,n_s_vec{3}*2,n_s_vec{4}};
% N_finite_elements_vec = 2:2:12;
N_finite_elements_vec = 2:7;
N_experiments = 1:4;
%% for storing results
f_obj_fesd  = {};
cpu_fesd  = {};
comp_resiudal_fesd = {};
error_temrinal_fesd = {};
%
f_obj_std  = {};
cpu_std  = {};
comp_resiudal_std = {};
error_temrinal_std = {};

%% FESD
for k = 1:length(N_experiments)
    settings.irk_scheme = irk_schemes{k};

    f_obj_fesd_k  = [];
    cpu_fesd_k  = [];
    comp_resiudal_fesd_k = [];
    error_temrinal_fesd_k = [];
    % std
    f_obj_std_k  = [];
    cpu_std_k  = [];
    comp_resiudal_std_k = [];
    error_temrinal_std_k = [];
    settings.use_fesd = 1;
    settings.mpcc_mode= mpcc_mode_fesd(k);
    for jj = n_s_vec{k}
        f_obj_temp = [];
        cpu_temp = [];
        comp_resiudal_temp = [];
        error_temrinal_temp = [];
        for ii = N_finite_elements_vec
            settings.n_s = jj;
            scenario.N_finite_elements = ii;
            scenario.scenario_name = [settings.irk_scheme '_n_s_' num2str(settings.n_s) '_N_FE_' num2str(scenario.N_finite_elements) '_FESD_' num2str(settings.use_fesd)];

            [results] = sliding_mode_ocp_experiment(scenario,model,settings);
            f_obj_temp  = [f_obj_temp,results.f_opt];
            cpu_temp  = [cpu_temp,results.cpu_time];
            comp_resiudal_temp = [comp_resiudal_temp, results.comp_residual];
            error_temrinal_temp = [error_temrinal_temp, results.error_terminal];
        end
        f_obj_fesd_k  = [f_obj_fesd_k;f_obj_temp];
        cpu_fesd_k  = [cpu_fesd_k ;cpu_temp];
        comp_resiudal_fesd_k = [comp_resiudal_fesd_k ;comp_resiudal_temp];
        error_temrinal_fesd_k = [error_temrinal_fesd_k;error_temrinal_temp];
    end

    % std experiments
    settings.use_fesd = 0;
    settings.mpcc_mode= mpcc_mode_std(k);
    for jj = n_s_vec{k}
        f_obj_temp = [];
        cpu_temp = [];
        comp_resiudal_temp = [];
        error_temrinal_temp = [];
        for ii = N_finite_elements_vec
            settings.n_s = jj;
            scenario.N_finite_elements = ii;
            scenario.scenario_name = [settings.irk_scheme '_n_s_' num2str(settings.n_s) '_N_FE_' num2str(scenario.N_finite_elements) '_FESD_' num2str(settings.use_fesd)];
            [results] = sliding_mode_ocp_experiment(scenario,model,settings);
            f_obj_temp  = [f_obj_temp,results.f_opt];
            cpu_temp  = [cpu_temp,results.cpu_time];
            comp_resiudal_temp = [comp_resiudal_temp,results.comp_residual];
            error_temrinal_temp = [error_temrinal_temp, results.error_terminal];
        end
        f_obj_std_k  = [f_obj_std_k;f_obj_temp];
        cpu_std_k  = [cpu_std_k ;cpu_temp];
        comp_resiudal_std_k = [comp_resiudal_std_k ;comp_resiudal_temp];
        error_temrinal_std_k = [error_temrinal_std_k;error_temrinal_temp];
    end
    % store all data
    f_obj_fesd  = [f_obj_fesd;f_obj_fesd_k  ];
    cpu_fesd  = [cpu_fesd ;cpu_fesd_k  ];
    comp_resiudal_fesd = [comp_resiudal_fesd;comp_resiudal_fesd_k];
    error_temrinal_fesd = [error_temrinal_fesd;error_temrinal_fesd_k];

    f_obj_std  = [f_obj_std;f_obj_std_k];
    cpu_std  = [cpu_std;cpu_std_k];
    comp_resiudal_std = [comp_resiudal_std;comp_resiudal_std_k];
    error_temrinal_std = [error_temrinal_std;error_temrinal_std_k];
end

%% Plot objective error
list_of_markers = {'-o','-d','-s','-x','-v','-^','-<','->','-*','-+','-|','-p','-h','-.'};
% f_star = f_obj_fesd_k(:,end);
f_star = f_obj_fesd{1}(end,end);
legend_str = {};
marker_counter = 1;
% Compare objectives
figure
for k = 1:length(N_experiments)
    error_obj_fesd = abs(f_obj_fesd{k}-f_star)+1e-16;
    subplot(121)
    for ii = 1:length(n_s_vec{k})
        %     semilogy(N_finite_elements_vec*N_stages,error_obj_fesd(ii,:),list_of_markers{ii},'linewidth',1.5);
        loglog(N_finite_elements_vec*N_stages,error_obj_fesd(ii,:),list_of_markers{ii},'linewidth',1.5);
        %     %     loglog(N_finite_elements_vec,error_obj_fesd(ii,:),list_of_markers{ii},'linewidth',1.5);
        hold on
        legend_str = [legend_str,[irk_schemes{k} '-' num2str(orders_vec{k}(ii))]];
        marker_counter  = marker_counter +1;
    end
    xlabel('$N_{\mathrm{FE}}$','interpreter','latex');
    ylabel('$|J_{\mathrm{FESD}}-J^*|$','interpreter','latex');
    grid on
    legend(legend_str,'Interpreter','latex','Location','best');
    ylim([1e-16 100])
    % xlim([2 N_finite_elements_vec(end)])
    error_obj_std= abs(f_obj_std{k}-f_star)+1e-16;
    % figure
    subplot(122)
    for ii = 1:length(n_s_vec{k})
        %     semilogy(N_finite_elements_vec*N_stages,error_obj_std(ii,:),list_of_markers{ii},'linewidth',1.5);
        loglog(N_finite_elements_vec*N_stages,error_obj_std(ii,:),list_of_markers{ii},'linewidth',1.5);
        hold on
    end
    xlabel('$N_{\mathrm{FE}}$','interpreter','latex');
    ylabel('$|J_{\mathrm{std}}-J^*|$','interpreter','latex');
    grid on
end
ylim([1e-16 100])
% xlim([2 N_finite_elements_vec(end)]);
saveas(gcf,'objective_error')

%% Plot Terminal state error (in one and two plot)
matlab_red = [0.8500 0.3250 0.0980];
matlab_blue  =[0 0.4470 0.7410];
legend_str = {};
marker_counter = 1;
if 0
    figure
    for k = 1:length(N_experiments)
        for ii = 1:length(n_s_vec{k})
            loglog(N_finite_elements_vec*N_stages,error_temrinal_fesd{k}(ii,:),list_of_markers{marker_counter},'linewidth',1.5);
            legend_str = [legend_str,[irk_schemes{k} '-' num2str(orders_vec{k}(ii)) '-FESD']];
            marker_counter  = marker_counter +1;
            hold on
            loglog(N_finite_elements_vec*N_stages,error_temrinal_std{k}(ii,:),list_of_markers{marker_counter},'linewidth',1.5);
            legend_str = [legend_str,[irk_schemes{k} '-' num2str(orders_vec{k}(ii)) '-Std']];
            grid on
            marker_counter  = marker_counter +1;
        end
        legend(legend_str,'Interpreter','latex','Location','best');
    end
    saveas(gcf,'terminal_error_1')
end
%%
list_of_markers = {'-o','-d','-s','-x','-v','-^','-<','->','-*','-+','-|','-p','-h','-.'};

figure
legend_str = {};
marker_counter = 1;
for k = 1:length(N_experiments)
    %     figure
    for ii = 1:length(n_s_vec{k})
        subplot(121)
        loglog(N_finite_elements_vec*N_stages,error_temrinal_fesd{k}(ii,:),list_of_markers{marker_counter},'linewidth',1.5);
        %         legend_str = [legend_str,[irk_schemes{k} '-' num2str(orders_vec(ii)) '-FESD']];
        %         marker_counter  = marker_counter +1;
        hold on
        grid on
        subplot(122)
        loglog(N_finite_elements_vec*N_stages,error_temrinal_std{k}(ii,:),list_of_markers{marker_counter},'linewidth',1.5);
        legend_str = [legend_str,[irk_schemes{k} ' ' num2str(orders_vec{k}(ii))]];
        grid on
        hold on
        marker_counter  = marker_counter +1;
    end
    subplot(121)
    legend(legend_str,'Interpreter','latex','Location','best');
    ylim([1e-16 1e1 ])
    subplot(122)
    ylim([1e-16 1e1 ])
end
xlabel('$N_{\mathrm{FE}}$','interpreter','latex');
ylabel('Terminal Error','interpreter','latex');
saveas(gcf,'terminal_error_2')

%%
list_of_markers = {'-o','-d','-s','-x','-v','-^','-<','->','-*','-+','-|','-p','-h','-.'};
for k = 1:length(N_experiments)
    figure
    legend_str = {};
    marker_counter = 1;
    for ii = 1:length(n_s_vec{k})
        subplot(121)
        loglog(N_finite_elements_vec*N_stages,error_temrinal_fesd{k}(ii,:),list_of_markers{marker_counter},'linewidth',1.5);
        hold on
        grid on
        subplot(122)
        loglog(N_finite_elements_vec*N_stages,error_temrinal_std{k}(ii,:),list_of_markers{marker_counter},'linewidth',1.5);
        legend_str = [legend_str,[irk_schemes{k} ' ' num2str(orders_vec{k}(ii))]];
        grid on
        hold on
        marker_counter  = marker_counter +1;
    end
    subplot(121)
    legend(legend_str,'Interpreter','latex','Location','best');
    ylim([1e-16 1e1 ])
    subplot(122)
    ylim([1e-16 1e1 ])
    xlabel('$N_{\mathrm{FE}}$','interpreter','latex');
    ylabel('Terminal Error','interpreter','latex');
    saveas(gcf,['terminal_error_' irk_schemes{k}])
end


%% Plot error as function of CPU time (pareto plot)
figure
list_of_markers = {'o','d','s','x','v','^','<','>','*','+','|','p','h','.'};

marker_counter = 1;
legend_str = {};
for k = 1:length(N_experiments)
    for ii = 1:length(n_s_vec{k})
        semilogy(cpu_fesd{k}(ii,:),error_temrinal_fesd{k}(ii,:),list_of_markers{marker_counter},'linewidth',1.5,'Color',matlab_red);
        legend_str = [legend_str,[irk_schemes{k} '-' num2str(orders_vec{k}(ii)) '-FESD']];
        %         marker_counter  = marker_counter +1;
        hold on
        semilogy(cpu_std{k}(ii,:),error_temrinal_std{k}(ii,:),list_of_markers{marker_counter},'linewidth',1.5,'Color',matlab_blue);
        legend_str = [legend_str,[irk_schemes{k} '-' num2str(orders_vec{k}(ii)) '-Std']];
        grid on
        marker_counter  = marker_counter+1;
    end
end
legend(legend_str,'Interpreter','latex','Location','best');
xlabel('CPU Time [s]','interpreter','latex');
ylabel('Terminal Error','interpreter','latex');
saveas(gcf,'pareto_fesd_vs_std')

%%
figure
list_of_markers = {'o','d','s','x','v','^','<','>','*','+','|','p','h','.'};
marker_counter = 1;
legend_str = {};
for k = 1:length(N_experiments)
    legend_str = [legend_str,[irk_schemes{k} '-FESD']];
    legend_str = [legend_str,[irk_schemes{k} '-Std']];
    semilogy(cpu_fesd{k}(:),error_temrinal_fesd{k}(:),list_of_markers{marker_counter},'linewidth',1.5,'Color',matlab_red);
    hold on
    semilogy(cpu_std{k}(:),error_temrinal_std{k}(:),list_of_markers{marker_counter},'linewidth',1.5,'Color',matlab_blue);

    marker_counter  = marker_counter+1;
    legend(legend_str,'Interpreter','latex','Location','best');
end
grid on
xlabel('CPU Time [s]','interpreter','latex');
ylabel('Terminal Error','interpreter','latex');
saveas(gcf,'pareto_fesd_vs_std_all_together')
%%
figure
list_of_markers = {'o','d','s','p','x','v','^','<','>','*'};
marker_counter = 1;
legend_str = {};
line_counter = 0;
include_lines_in_legend = [];
plot_handles = {};
for k = 1:length(N_experiments)
    for ii = 1:length(n_s_vec{k})
        line_counter = line_counter+1;

        plot_handles{line_counter} = semilogy(cpu_fesd{k}(ii,:),error_temrinal_fesd{k}(ii,:),list_of_markers{marker_counter},'linewidth',1.5,'Color',matlab_red,'MarkerSize',2+(ii*2+1));
        if ii == 1
            legend_str = [legend_str,[irk_schemes{k} '-FESD']];
            include_lines_in_legend = [include_lines_in_legend,line_counter];
        else
            legend_str = [legend_str,'x'];
        end
        hold on
        line_counter = line_counter+1;
        plot_handles{line_counter} =  semilogy(cpu_std{k}(ii,:),error_temrinal_std{k}(ii,:),list_of_markers{marker_counter},'linewidth',1.5,'Color',matlab_blue,'MarkerSize',2+(ii*2+1));
        if ii == 1
            legend_str = [legend_str,[irk_schemes{k} '-Std']];
            include_lines_in_legend = [include_lines_in_legend,line_counter];
        else
            legend_str = [legend_str,'x'];
        end
    end
    marker_counter  = marker_counter+1;
end
grid on
% legend(legend_str,'Interpreter','latex','Location','best')

hleg = legend([plot_handles{include_lines_in_legend}],'Interpreter','latex');
hleg.String = legend_str(include_lines_in_legend);
hleg.Location = 'southeast';

xlabel('CPU Time [s]','interpreter','latex');
ylabel('Terminal Error','interpreter','latex');
ylim([1e-12 10])
xlim([0 30])
saveas(gcf,'pareto_fesd_vs_std_cpu_lin')
xlim([0 250])
set(gca, 'XScale', 'log')
saveas(gcf,'pareto_fesd_vs_std_cpu_log')
%% pareto plot as function of number of stage points
figure
list_of_markers = {'o','d','s','p','x','v','^','<','>','*'};
marker_counter = 1;
legend_str = {};
line_counter = 0;
include_lines_in_legend = [];
plot_handles = {};
for k = 1:length(N_experiments)
    for ii = 1:length(n_s_vec{k})
        line_counter = line_counter+1;
        plot_handles{line_counter} = semilogy(n_s_vec{k}(ii)*N_finite_elements_vec*N_stages,error_temrinal_fesd{k}(ii,:),list_of_markers{marker_counter},'linewidth',1.5,'Color',matlab_red,'MarkerSize',3+(ii*2+1));
        if ii == 1
            legend_str = [legend_str,[irk_schemes{k} '-FESD']];
            include_lines_in_legend = [include_lines_in_legend,line_counter];
        else
            legend_str = [legend_str,'x'];
        end
        hold on
        line_counter = line_counter+1;
        plot_handles{line_counter} =semilogy(n_s_vec{k}(ii)*N_finite_elements_vec*N_stages,error_temrinal_std{k}(ii,:),list_of_markers{marker_counter},'linewidth',1.5,'Color',matlab_blue,'MarkerSize',3+(ii*2+1));
        if ii == 1
            legend_str = [legend_str,[irk_schemes{k} '-Std']];
            include_lines_in_legend = [include_lines_in_legend,line_counter];
        else
            legend_str = [legend_str,'x'];
        end
    end
    marker_counter  = marker_counter+1;
end
grid on
% legend(legend_str{1:length(N_experiments)},'Interpreter','latex','Location','best')

hleg = legend([plot_handles{include_lines_in_legend}],'Interpreter','latex');
hleg.String = legend_str(include_lines_in_legend);
hleg.Location = 'best';

xlabel('Total Number of Stage Points','interpreter','latex');
ylabel('Terminal Error','interpreter','latex');
ylim([1e-12 10])
saveas(gcf,'pareto_fesd_vs_std_stage_pt_lin')
set(gca, 'XScale', 'log')
saveas(gcf,'pareto_fesd_vs_std_stage_pt_log')

%%
% results std
results_benchmark.f_obj_std  = f_obj_std ;
results_benchmark.cpu_std  = cpu_std ;
results_benchmark.list_of_markers = list_of_markers;
results_benchmark.comp_resiudal_std  = comp_resiudal_std;
results_benchmark.error_temrinal_std  = error_temrinal_std;
% results fesd
results_benchmark.f_obj_fesd  = f_obj_fesd;
results_benchmark.cpu_fesd  = cpu_fesd ;
results_benchmark.comp_resiudal_fesd  = comp_resiudal_fesd;
results_benchmark.error_temrinal_fesd  = error_temrinal_fesd;
% scnearios
results_benchmark.irk_schemes = irk_schemes;
results_benchmark.N_experiments = N_experiments;
results_benchmark.N_finite_elements_vec  = N_finite_elements_vec ;
results_benchmark.n_s_vec  = n_s_vec ;

if scenario.save_results
    save(['results/results_benchmark_full.mat'],'results_benchmark')
    %     save([scenario_name '.mat'],'results')
end