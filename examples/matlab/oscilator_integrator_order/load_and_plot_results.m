% clear all
% close all
%%
all_scenarios = {'irk_radau_iia_fesd_differential'};
load([all_scenarios{1} '.mat'])
unfold_struct(results,'caller')

%% Error plots
list_of_markers = {'-o','-d','-s','-x','-v','-^','-<','->','-*'};
% error as function of step size
figure
for ii = 1:length(n_s_vec)
    loglog(nominal_h_all_experiments(ii,:),complementarity_all_experiments(ii,:)+1e-18,list_of_markers{ii},'linewidth',1.5);
    hold on
end
xlabel('$h$','interpreter','latex');
ylabel('comp residual','interpreter','latex');
grid on
legend(legend_str,'interpreter','latex');
% as function of num of stage points
figure
for ii = 1:length(n_s_vec)
    loglog(M_true_all_experiment(ii,:),errors_all_experiments(ii,:),list_of_markers{ii},'linewidth',1.5);
    hold on
end
xlabel('$M$','interpreter','latex');
ylabel('$E(2)$','interpreter','latex');
grid on
ylim([1e-15 100])
legend(legend_str,'interpreter','latex');
%% complementarity residual
% as function of step size
figure
for ii = 1:length(n_s_vec)
    loglog(nominal_h_all_experiments(ii,:),errors_all_experiments(ii,:),list_of_markers{ii},'linewidth',1.5);
    hold on
    xlabel('$h$','interpreter','latex');
    ylabel('$E(2)$','interpreter','latex');
    grid on
end
ylim([1e-15 100])
legend(legend_str,'interpreter','latex');

figure
loglog(M_true_all_experiment,complementarity_all_experiments+1e-18,list_of_markers{ii},'linewidth',1.5);
hold on
xlabel('$M$','interpreter','latex');
ylabel('comp residual','interpreter','latex');
grid on
legend(legend_str,'interpreter','latex');

%% Practial slopes
practical_slopes_all_experiments = [];
for ii = 1:length(n_s_vec)
    practical_slopes_all_experiments = [practical_slopes_all_experiments;diff(log(errors_all_experiments(ii,:)))./diff(log(nominal_h_all_experiments(ii,:)))];
end

figure
bar([practical_slopes_all_experiments,mean(practical_slopes_all_experiments,2)])
xticks(1:length(n_s_vec))
xticklabels(legend_str);
grid on
