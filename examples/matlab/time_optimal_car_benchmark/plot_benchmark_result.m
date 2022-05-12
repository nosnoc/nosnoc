% load results
load('output_fesd.mat')
markers = {'p','s','d','h'};
N_stages_vec = output_fesd.N_stages_vec;
run_experiments = output_fesd.run_experiments;
N_trails = output_fesd.N_trails;
% experiment_names = output_fesd.experiment_names;
legend_str = output_fesd.legend_str;
legend_str_temp = legend_str;
legend_str = {};
cpu_time = [];
error = [];
T_opt = [];
% read all data
for ii = 1:length(run_experiments)
        if run_experiments(ii) == 1
            load(['output_' experiment_names{ii} '.mat'])
            eval(['unfold_struct(output_' experiment_names{ii} ',''caller'')']);
            legend_str  = [legend_str,legend_str_temp{ii}];
            eval(['error = [error;error_' experiment_names{ii} '];']);
            eval(['cpu_time = [cpu_time;cpu_time_' experiment_names{ii} '];']);
            eval(['T_opt  = [T_opt ;T_opt_' experiment_names{ii} '];']);
        end
end

%% plot cpu time linear scale
figure
for ii = 1:sum(run_experiments)
    plot(N_stages_vec,cpu_time(ii,:),[markers{ii} '-'],'LineWidth',1.5);
    hold on
end
xlabel('$N_{\mathrm{stages}}$','Interpreter','latex');
ylabel('CPU Time [s]','Interpreter','latex');
legend(legend_str,'Interpreter','latex','Location','best');
grid on
%% plot cpu time semilogy scale
figure
for ii = 1:sum(run_experiments)
    semilogy(N_stages_vec,cpu_time(ii,:),[markers{ii} '-'],'LineWidth',1.5);
    hold on
end
xlabel('$N_{\mathrm{stages}}$','Interpreter','latex');
ylabel('CPU Time [s]','Interpreter','latex');
legend(legend_str,'Interpreter','latex','Location','best');
grid on


figure
for ii = 1:sum(run_experiments)
    plot(N_stages_vec,T_opt(ii,:),[markers{ii} '-'],'LineWidth',1.5);
    hold on
end
xlabel('$N_{\mathrm{stages}}$','Interpreter','latex');
ylabel('$T^{\star}$','Interpreter','latex');
legend(legend_str,'Interpreter','latex','Location','best');
grid on
%% log log cput vs error
figure
for ii = 1:sum(run_experiments)
    loglog(cpu_time(ii,:),error(ii,:),[markers{ii}],'LineWidth',1.5);
    hold on
end
xlabel('CPU Time [s]','Interpreter','latex');
ylabel('$E(T^{\star})$','Interpreter','latex');
ylim([1e-15 5e2])
grid on
legend(legend_str,'Interpreter','latex','Location','best');
%% simlogy cput vs error
figure
for ii = 1:sum(run_experiments)
    semilogy(cpu_time(ii,:),error(ii,:),[markers{ii}],'LineWidth',1.5);
    hold on
end
xlabel('CPU Time [s]','Interpreter','latex');
ylabel('$E(T^{\star})$','Interpreter','latex');
ylim([1e-15 5e2])
grid on
legend(legend_str,'Interpreter','latex','Location','best');
%%
figure
for ii = 1:sum(run_experiments)
    loglog(N_stages_vec,error(ii,:),[markers{ii} '-'],'LineWidth',1.5);
    hold on
end
xlabel('$N_{\mathrm{stages}}$','Interpreter','latex');
ylabel('$E(T^{\star})$','Interpreter','latex');
ylim([1e-15 5e2])
grid on
legend(legend_str,'Interpreter','latex','Location','best');