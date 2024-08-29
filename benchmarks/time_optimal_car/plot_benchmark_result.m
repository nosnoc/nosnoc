%% Meta data
markers = {'p','s','d','h'};
N_stages_vec = [7:3:13];
N_stages_vec = [10:5:80];
N_trails = 1;
minlp_time_limit = 600;
experiment_names = {'fesd','std','gurobi','bonmin'};
legend_str = {'NOSNOC-FESD','NOSNOC-Std','Gurobi','Bonmin'};
run_fesd = 1;
run_std = 1;
run_gurobi = 0;
run_bonmin = 0;
run_experiments = [run_fesd run_std run_gurobi run_bonmin];
legend_str_temp = legend_str;
%%
% load(['output_gurobi.mat'])
%% Load results

legend_str = {};
cpu_time = [];
error = [];
T_opt = [];
% read all data
for ii = 1:length(run_experiments)
        if run_experiments(ii) == 1
            load(['output_' experiment_names{ii} '.mat'])
            eval(['unfold_struct(output_' experiment_names{ii} ',''caller'')']);
            % 

            eval(['error_temp = error_' experiment_names{ii} ';']);
            eval(['cpu_time_temp = cpu_time_' experiment_names{ii} ';']);
            eval(['T_opt_temp = T_opt_' experiment_names{ii} ';']);
            error = [error;[error_temp,nan*ones(1,length(N_stages_vec)-length(error_temp))]];
            cpu_time = [cpu_time;[cpu_time_temp,nan*ones(1,length(N_stages_vec)-length(cpu_time_temp))]];
            T_opt = [T_opt;[T_opt_temp,nan*ones(1,length(N_stages_vec)-length(T_opt_temp))]];
            
%             eval(['error = [error;error_' experiment_names{ii} '];']);
% %             eval(['cpu_time = [cpu_time;cpu_time_' experiment_names{ii} '];']);
%             eval(['T_opt  = [T_opt ;T_opt_' experiment_names{ii} '];']);
            %
            legend_str  = [legend_str,legend_str_temp{ii}];
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
subplot(121)
for ii = 1:sum(run_experiments)
    semilogy(N_stages_vec,cpu_time(ii,:),[markers{ii} '-'],'LineWidth',1.5);
    hold on
end
xlabel('$N_{\mathrm{stages}}$','Interpreter','latex');
ylabel('CPU Time [s]','Interpreter','latex');
legend(legend_str,'Interpreter','latex','Location','best');
grid on

subplot(122)
for ii = 1:sum(run_experiments)
    plot(N_stages_vec,T_opt(ii,:),[markers{ii} '-'],'LineWidth',1.5);
    hold on
end
xlabel('$N_{\mathrm{stages}}$','Interpreter','latex');
ylabel('$T^{\star}$','Interpreter','latex');
legend(legend_str,'Interpreter','latex','Location','best');
grid on

%% cpu and error
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultlegendinterpreter','latex')
set(groot,'DefaultTextarrowshapeInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');
figure
subplot(121)
for ii = 1:sum(run_experiments)
    semilogy(N_stages_vec,cpu_time(ii,:),[markers{ii} '-'],'LineWidth',1.5);
    hold on
end
xlabel('$N_{\mathrm{stages}}$','Interpreter','latex');
ylabel('CPU Time [s]','Interpreter','latex');
% legend(legend_str,'Interpreter','latex','Location','best');
grid on

subplot(122)
for ii = 1:sum(run_experiments)
    loglog(cpu_time(ii,:),error(ii,:),[markers{ii}],'LineWidth',1.5);
    hold on
end
xlabel('CPU Time [s]','Interpreter','latex');
ylabel('$E(T)$','Interpreter','latex');
ylim([1e-15 5e2])
grid on
legend(legend_str,'Interpreter','latex','Location','best');


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
