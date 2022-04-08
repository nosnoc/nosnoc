% clear all
close all
%%
all_scenarios_fesd = {'irk_radau_iia_fesd_differential','irk_gauss_legendre_fesd_differential',...
                    'irk_radau_iia_fesd_integral','irk_gauss_legendre_fesd_integral',...
                        'irk_radau_ia_fesd_differential','irk_radau_i_fesd_differential', ...
                       'irk_Lobatto_iii_fesd_differential','irk_Lobatto_iiia_fesd_differential',...
                       'irk_Lobatto_iiib_fesd_differential','irk_Lobatto_iiic_fesd_differential','erk_fesd_differential'};
all_scenarios_std =  {'irk_radau_iia_std_differential','irk_gauss_legendre_std_differential',...
                    'irk_radau_iia_std_integral','irk_gauss_legendre_std_integral',...
                        'irk_radau_ia_std_differential','irk_radau_i_std_differential', ...
                       'irk_Lobatto_iii_std_differential','irk_Lobatto_iiia_std_differential',...
                       'irk_Lobatto_iiib_std_differential','irk_Lobatto_iiic_std_differential','erk_std_differential'};

scenario_num = 10;
partical_row = 1;
list_of_markers = {'-o','-d','-s','-x','-v','-^','-<','->','-*'};
fprintf(['Scenarios: ' all_scenarios_fesd{scenario_num} '.\n']);


%% Joint plot error


figure
subplot(121)
load([all_scenarios_std{scenario_num} '.mat'])
unfold_struct(results,'caller')
consider_methods = length(results.n_s_vec);
% error as function of step size
for ii = 1:length(n_s_vec)
    loglog(nominal_h_all_experiments(ii,1:consider_methods),errors_all_experiments(ii,1:consider_methods),list_of_markers{ii},'linewidth',1.5);
    hold on
end
xlabel('$h$','interpreter','latex');
ylabel('$E(T)$','interpreter','latex');
grid on
legend(legend_str{1:consider_methods},'interpreter','latex','Location','southwest');
ylim([1e-15 100])
subplot(122)
load([all_scenarios_fesd{scenario_num} '.mat'])
unfold_struct(results,'caller')
for ii = 1:length(n_s_vec)
    loglog(nominal_h_all_experiments(ii,1:consider_methods),errors_all_experiments(ii,1:consider_methods),list_of_markers{ii},'linewidth',1.5);
    hold on
end
xlabel('$h$','interpreter','latex');
grid on
ylim([1e-15 100])
%% Complementarity residual as sanity check

figure
subplot(121)
load([all_scenarios_std{scenario_num} '.mat'])
unfold_struct(results,'caller')
consider_methods = length(results.n_s_vec);
% error as function of step size
for ii = 1:length(n_s_vec)
    loglog(nominal_h_all_experiments(ii,1:consider_methods),complementarity_all_experiments(ii,1:consider_methods),list_of_markers{ii},'linewidth',1.5);
    hold on
end
xlabel('$h$','interpreter','latex');
ylabel('complementarity residual','interpreter','latex');
grid on
legend(legend_str{1:consider_methods},'interpreter','latex','Location','northwest');
ylim([1e-20 1])
subplot(122)
load([all_scenarios_fesd{scenario_num} '.mat'])
unfold_struct(results,'caller')
for ii = 1:length(n_s_vec)
    loglog(nominal_h_all_experiments(ii,1:consider_methods),complementarity_all_experiments(ii,1:consider_methods),list_of_markers{ii},'linewidth',1.5);
    hold on
end
xlabel('$h$','interpreter','latex');
grid on
ylim([1e-20 1])

%% Practial slopes
figure
subplot(121)
load([all_scenarios_std{scenario_num} '.mat'])
unfold_struct(results,'caller')
% bar([practical_slopes_all_experiments,mean(practical_slopes_all_experiments,2)])
% bar(max(practical_slopes_all_experiments(1:consider_methods),1))
bar([practical_slopes_all_experiments(1:consider_methods,partical_row)])
bar(mean(practical_slopes_all_experiments,2))
xticks(1:consider_methods)
xticklabels(legend_str);
ylabel('Estimate of integration order','Interpreter','latex');
grid on
ylim([0 10])
subplot(122)
load([all_scenarios_fesd{scenario_num} '.mat'])
unfold_struct(results,'caller')
% bar([practical_slopes_all_experiments,mean(practical_slopes_all_experiments,2)])
bar([practical_slopes_all_experiments(1:consider_methods,partical_row)])
% bar((practical_slopes_all_experiments(1:consider_methods),1))
% bar(mean(practical_slopes_all_experiments,2))
xticks(1:consider_methods)
xticklabels(legend_str);
grid on
ylim([0 10])

% best slope
figure
subplot(121)
load([all_scenarios_std{scenario_num} '.mat'])
unfold_struct(results,'caller')
bar(max(practical_slopes_all_experiments(1:consider_methods),1))

xticks(1:consider_methods)
xticklabels(legend_str);
ylabel('Estimate of integration order','Interpreter','latex');
grid on
ylim([0 10])
subplot(122)
load([all_scenarios_fesd{scenario_num} '.mat'])
unfold_struct(results,'caller')
bar(max(practical_slopes_all_experiments(1:consider_methods),1))
xticks(1:consider_methods)
xticklabels(legend_str);
grid on
ylim([0 10])

%% as function of num of stage points
% figure
% for ii = 1:length(n_s_vec)
%     loglog(M_true_all_experiment(ii,1:consider_methods),errors_all_experiments(ii,1:consider_methods),list_of_markers{ii},'linewidth',1.5);
%     hold on
% end
% xlabel('$M$','interpreter','latex');
% ylabel('$E(2)$','interpreter','latex');
% grid on
% ylim([1e-15 100])
% legend(legend_str,'interpreter','latex');
%% complementarity residual
% % as function of step size
% figure
% for ii = 1:length(n_s_vec)
%     loglog(nominal_h_all_experiments(ii,1:consider_methods),errors_all_experiments(ii,1:consider_methods),list_of_markers{ii},'linewidth',1.5);
%     hold on
%     xlabel('$h$','interpreter','latex');
%     ylabel('$E(2)$','interpreter','latex');
%     grid on
% end
% ylim([1e-15 100])
% legend(legend_str,'interpreter','latex');
% 
% figure
% loglog(M_true_all_experiment,complementarity_all_experiments+1e-18,list_of_markers{ii},'linewidth',1.5);
% hold on
% xlabel('$M$','interpreter','latex');
% ylabel('comp residual','interpreter','latex');
% grid on
% legend(legend_str,'interpreter','latex');