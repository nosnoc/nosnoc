%% ===================== SCRIPT 2: LOAD + COMPARE DATASETS =====================
clear; clc; close all;

latexify_plot()
nice_plot_colors;
linewidth = 1.5;

data_names = {'ideal_mpc','rti_mpc','as_mpc_qpcc','as_mpc_qp'};   % first entry is the reference solution
data_names =  {'Ideal-MPC','HyRTI','HyAS-RTI-QPCC','HyAS-RTI-QP','Smoothing','HyRTI-fast'};
data_names = {'Ideal-MPC','HyRTI','Smoothing','HyAS-RTI-QPCC'};
data_names = {'HyAS-RTI-QPCC'};
data_names = {'Full-MPCC-MPC','HyRTI', 'HyAS-RTI-QPCC', 'HyAS-RTI-QPCC-FULL', 'Smoothing'};

% line styles (cycled)
ls = {'-','--',':','-.'};

% load all
D = cell(size(data_names));
for k = 1:numel(data_names)
    S = load([data_names{k} '.mat'], 'data');
    D{k} = S.data;
end
Dref = D{1};

leg = cell(size(D));
for k = 1:numel(D)
    if isfield(D{k},'scenario')
        leg{k} = D{k}.scenario;
    else
        leg{k} = D{k}.name;
    end
end
%% Plot controls
ax1 = subplot(2,1,1);
hold on
ax2 = subplot(2,1,2);
hold on
h_1 = gobjects(numel(D),1);
h_2 = gobjects(numel(D),1);

for k = 1:numel(D)
    tk = D{k}.t_u;
    uk = D{k}.u;
    h_1(k) = stairs(ax1, tk, [uk(1,:), uk(1,end)]', 'LineWidth', linewidth, ...
        'LineStyle', ls{mod(k-1,numel(ls))+1});
    h_2(k) = stairs(ax2, tk, [uk(2,:), uk(2,end)]', 'LineWidth', linewidth, ...
        'LineStyle', ls{mod(k-1,numel(ls))+1});
end
xlabel('$t$')
ylabel('$u$')
grid on
legend(h_1, leg, 'Interpreter','none', 'Location','best');
legend(h_2, leg, 'Interpreter','none', 'Location','best');

%% Plot Computation times
figure;

subplot(121); hold on
hp = gobjects(numel(D),1);
for k = 1:numel(D)
    hp(k) = stairs(D{k}.t_u(1:end-1), D{k}.preparation_times, ...
        'LineWidth', linewidth, 'LineStyle', ls{mod(k-1,numel(ls))+1});
end
yline(D{1}.DT,'-','LineWidth',linewidth,'color',0.5*ones(3,1))
xlabel('$t$'); ylabel('Preparation time [s]'); grid on
set(gca,'YScale','log')
legend(hp, leg, 'Interpreter','none', 'Location','best');

subplot(122); hold on
hf = gobjects(numel(D),1);
for k = 1:numel(D)
    hf(k) = stairs(D{k}.t_u(1:end-1), D{k}.feedback_times, ...
        'LineWidth', linewidth, 'LineStyle', ls{mod(k-1,numel(ls))+1});
end
yline(D{1}.DT,'-','LineWidth',linewidth,'color',0.5*ones(3,1))
xlabel('$t$'); ylabel('Feedback time [s]'); grid on
set(gca,'YScale','log')
legend(hf, leg, 'Interpreter','none', 'Location','best');

%% State trajectories
figure;

ls = {'-','--',':','-.'}; % if not defined above

% common settings (assume comparable)
u_max = D{1}.u_max;

for k = 1:numel(D)
    stylek = ls{mod(k-1,numel(ls))+1};

    % use high-res for states, MPC grid for u
    t   = D{k}.t_high_res;
    x   = D{k}.x_high_res;
    t_u = D{k}.t_u;
    u   = D{k}.u;

    DT = D{k}.DT;
    theta_ref = D{k}.x_refs(3,:);
    
    % detect zero crossings of v = x(3,:)
    v = x(6,:);
    
    % --- q/theta
    subplot(411); hold on
    plot(t,x(1,:),'LineWidth',linewidth,'LineStyle',stylek,'DisplayName',[leg{k} ' $q_1$'])
    plot(t,x(2,:),'LineWidth',linewidth,'LineStyle',stylek,'DisplayName',[leg{k} ' $q_2$'])
    xlabel("$t$"); ylabel("$q$")
    % --- v/omega + switch markers
    subplot(412); hold on
    plot(t,x(3,:),'LineWidth',linewidth,'LineStyle',stylek,'DisplayName',[leg{k} ' $\theta$'])
    % reference (draw once, not per dataset)
    if k == 1
        stairs(t_u(1:end-1), theta_ref,'k--','LineWidth',linewidth,'HandleVisibility','off')
        xlabel("$t$"); ylabel("$\theta$")
    end
    % --- omega
    subplot(413); hold on
    plot(t,x(6,:),'LineWidth',linewidth,'LineStyle',stylek,'DisplayName',[leg{k} ' $\omega$'])
    xlabel("$t$"); ylabel("$\omega$")
    % --- omega
    subplot(414); hold on
    plot(t,x(4,:),'LineWidth',linewidth,'LineStyle',stylek,'DisplayName',[leg{k} ' $\dot{q}_1$'])
    plot(t,x(5,:),'LineWidth',linewidth,'LineStyle',stylek,'DisplayName',[leg{k} ' $\dot{q}_2$'])
    xlabel("$t$"); ylabel("$\dot q$")
end

subplot(411); legend('Location','southeast');
subplot(412); legend('Location','southeast');
subplot(413); legend('Location','best');
subplot(414); legend('Location','best');
