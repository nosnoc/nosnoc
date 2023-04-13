% sliding_moode_ocp_plots
if 0
N_stages = 6;
mpcc_mode_fesd = [5 5 3 3];
mpcc_mode_std = [3 3 3 3];
irk_schemes = {IRKSchemes.RADAU_IIA,'Lobatto-IIIC','Gauss-Legendre','Explicit-RK'};
n_s_vec = {[1:4],[2:4],[1:4],[1 4]};
orders_vec = {n_s_vec{1}*2-1,n_s_vec{2}*2-2,n_s_vec{3}*2,n_s_vec{4}};
load(['results\results_benchmark_full.mat'])
unfold_struct(results_benchmark,'base')
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
xlabel('CPU time [s]','interpreter','latex');
ylabel('Terminal constraint satisfaction','interpreter','latex');

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
xlabel('CPU time [s]','interpreter','latex');
ylabel('Terminal constraint satisfaction','interpreter','latex');

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

        plot_handles{line_counter} = semilogy(cpu_fesd{k}(ii,:),error_temrinal_fesd{k}(ii,:),list_of_markers{marker_counter},'linewidth',1.5,'Color',matlab_red,'MarkerSize',2+(ii*2));
        if ii == 1
            legend_str = [legend_str,[irk_schemes{k} '-FESD']];
            include_lines_in_legend = [include_lines_in_legend,line_counter];
        else
            legend_str = [legend_str,'x'];
        end
        hold on
        line_counter = line_counter+1;
        plot_handles{line_counter} =  semilogy(cpu_std{k}(ii,:),error_temrinal_std{k}(ii,:),list_of_markers{marker_counter},'linewidth',1.5,'Color',matlab_blue,'MarkerSize',2+(ii*2));
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

xlabel('CPU time [s]','interpreter','latex');
ylabel('Terminal constraint satisfaction','interpreter','latex');

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
        plot_handles{line_counter} = semilogy(n_s_vec{k}(ii)*N_finite_elements_vec*N_stages,error_temrinal_fesd{k}(ii,:),list_of_markers{marker_counter},'linewidth',1.5,'Color',matlab_red,'MarkerSize',2+(ii*2));
        if ii == 1
            legend_str = [legend_str,[irk_schemes{k} '-FESD']];
            include_lines_in_legend = [include_lines_in_legend,line_counter];
        else
            legend_str = [legend_str,'x'];
        end
        hold on
        line_counter = line_counter+1;
        plot_handles{line_counter} =semilogy(n_s_vec{k}(ii)*N_finite_elements_vec*N_stages,error_temrinal_std{k}(ii,:),list_of_markers{marker_counter},'linewidth',1.5,'Color',matlab_blue,'MarkerSize',2+(ii*2));
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

xlabel('Total number of stages','interpreter','latex');
ylabel('Terminal constraint satisfaction','interpreter','latex');
ylim([1e-12 10])
saveas(gcf,'pareto_fesd_vs_std_stage_pt_lin')
set(gca, 'XScale', 'log')
saveas(gcf,'pareto_fesd_vs_std_stage_pt_log')