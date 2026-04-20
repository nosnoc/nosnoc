%% ===================== LOAD + COMPARE DATASETS =====================
clear; clc; close all;

latexify_plot()
nice_plot_colors;
linewidth = 1.5;

data_names = {'ideal_mpc','rti_mpc','as_mpc_qpcc','as_mpc_qp'};   % first entry is the reference solution
data_names =  {'Ideal-MPC','HyRTI','HyAS-RTI-QPCC','HyAS-RTI-QP','Smoothing','HyRTI-fast'};
data_names = {'Ideal-MPC','HyRTI', 'HyAS-RTI-QPCC'};
data_names = {'HyRTI', 'HyAS-RTI-QPCC','HyRTI-IPOPT', 'HyAS-RTI-QPCC-IPOPT',};


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
%% 1) Plot controls u(t) for each dataset (overlay)
figure; hold on
h = gobjects(numel(D),1);

for k = 1:numel(D)
    tk = D{k}.t_u;
    uk = D{k}.u;
    h(k) = stairs(tk, [uk, uk(end)], 'LineWidth', linewidth, ...
        'LineStyle', ls{mod(k-1,numel(ls))+1});
end
xlabel('$t$')
ylabel('$u$')
grid on
legend(h, leg, 'Interpreter','none', 'Location','best');

%% 3) semilogy of distance to piecewise setpoint for each dataset (using HIGH-RES state)
figure; hold on
h3 = gobjects(numel(D),1);

for k = 1:numel(D)
    tk = D{k}.t_high_res;
    xk = D{k}.x_high_res;

    t_switch = D{k}.t_new_set_point;

    xref = repmat(D{k}.x_ref1, 1, numel(tk));
    idx2 = tk > t_switch;
    xref(:,idx2) = repmat(D{k}.x_ref2, 1, nnz(idx2));

    dist = vecnorm(xk - xref, 2, 1);

    h3(k) = plot(tk, dist, 'LineWidth', linewidth, ...
        'LineStyle', ls{mod(k-1,numel(ls))+1});
end

xlabel('$t$')
ylabel('$\|x_{\mathrm{hi}} - x_{\mathrm{ref}}(t)\|_2$')
grid on
set(gca,'YScale','log')
legend(h3, leg, 'Interpreter','none', 'Location','best');

%% 6) value functions
figure; hold on
h3 = gobjects(numel(D),1);

for k = 1:numel(D)
    tk = D{k}.t_u;
    V = D{k}.f_opt;
    h3(k) = plot(tk,[V,nan], 'LineWidth', linewidth, 'LineStyle', ls{mod(k-1,numel(ls))+1});
end
set(gca,'Yscale','log');
ylim([1e-2 1e4])

xlabel('$t$')
ylabel('$V(x(t))$')
grid on
set(gca,'YScale','log')
legend(h3, leg, 'Interpreter','none', 'Location','best');



%% 3a) semilogy of distance to piecewise setpoint for each dataset (only position)
figure; hold on
h3 = gobjects(numel(D),1);
idx_state = 1;

for k = 1:numel(D)
    tk = D{k}.t_high_res;
    xk = D{k}.x_high_res;

    t_switch = D{k}.t_new_set_point;

    xref = repmat(D{k}.x_ref1, 1, numel(tk));
    idx2 = tk > t_switch;
    xref(:,idx2) = repmat(D{k}.x_ref2, 1, nnz(idx2));

    dist = vecnorm(xk(idx_state,:) - xref(idx_state,:), 2, 1);

    h3(k) = plot(tk, dist, 'LineWidth', linewidth, ...
        'LineStyle', ls{mod(k-1,numel(ls))+1});
end

xlabel('$t$')
ylabel('$\|x_{\mathrm{hi}} - x_{\mathrm{ref}}(t)\|_2$')
grid on
set(gca,'YScale','log')
legend(h3, leg, 'Interpreter','none', 'Location','best');

%% 4) Computation times: preparation_time and feedback_time (stairs)
figure;

subplot(211); hold on
hp = gobjects(numel(D),1);
for k = 1:numel(D)
    hp(k) = stairs(D{k}.t_u(1:end-1), D{k}.preparation_times, ...
        'LineWidth', linewidth, 'LineStyle', ls{mod(k-1,numel(ls))+1});
end
yline(D{1}.DT,'-','LineWidth',linewidth,'color',0.5*ones(3,1))
xlabel('$t$'); ylabel('Preparation time [s]'); grid on
set(gca,'YScale','log')
legend(hp, leg, 'Interpreter','none', 'Location','best');

subplot(212); hold on
hf = gobjects(numel(D),1);
for k = 1:numel(D)
    hf(k) = stairs(D{k}.t_u(1:end-1), D{k}.feedback_times, ...
        'LineWidth', linewidth, 'LineStyle', ls{mod(k-1,numel(ls))+1});
end
yline(D{1}.DT,'-','LineWidth',linewidth,'color',0.5*ones(3,1))
xlabel('$t$'); ylabel('Feedback time [s]'); grid on
set(gca,'YScale','log')
legend(hf, leg, 'Interpreter','none', 'Location','best');

%% 5) State trajectories (all datasets) in your updated style
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
    t_switch_set_point = D{k}.t_new_set_point;
    x_ref1 = D{k}.x_ref1;
    x_ref2 = D{k}.x_ref2;

    % detect zero crossings of v = x(3,:)
    v = x(3,:);
    eps_v = 1e-1;
    state = zeros(size(v));
    state(v >  eps_v) =  1;
    state(v < -eps_v) = -1;
    s = state;
    last = 0;
    for i = 1:numel(s)
        if s(i) ~= 0, last = s(i); else, s(i) = last; end
    end
    idx_cross = find(s(1:end-1).*s(2:end) < 0);
    t_cross = t(idx_cross);

    % split time for references
    i_sw = find(t <= t_switch_set_point, 1, 'last');
    if isempty(i_sw), i_sw = 1; end
    t1 = t(1:i_sw);
    t2 = t(i_sw:end);

    % --- q/theta
    subplot(311); hold on
    plot(t,x(1,:),'LineWidth',linewidth,'LineStyle',stylek,'DisplayName',[leg{k} ' $q$'])
    plot(t,x(2,:),'LineWidth',linewidth,'LineStyle',stylek,'DisplayName',[leg{k} ' $\theta$'])

    % reference (draw once, not per dataset)
    if k == 1
        plot(t1, x_ref1(1)*ones(size(t1)),'k--','LineWidth',linewidth,'HandleVisibility','off')
        plot(t1, x_ref1(2)*ones(size(t1)),'k--','LineWidth',linewidth,'HandleVisibility','off')
        plot(t2, x_ref2(1)*ones(size(t2)),'k--','LineWidth',linewidth,'HandleVisibility','off')
        plot(t2, x_ref2(2)*ones(size(t2)),'k--','LineWidth',linewidth,'HandleVisibility','off')
        xlabel("$t$"); ylabel("$q,\, \theta$")
    end

    % --- v/omega + switch markers
    subplot(312); hold on
    plot(t,x(3,:),'LineWidth',linewidth,'LineStyle',stylek,'DisplayName',[leg{k} ' $v$'])
    plot(t,x(4,:),'LineWidth',linewidth,'LineStyle',stylek,'DisplayName',[leg{k} ' $\omega$'])
    for j = 1:numel(t_cross)
        xline(t_cross(j),'k:','LineWidth',1.2,'HandleVisibility','off')
    end
    if k == 1
        xlabel("$t$"); ylabel("$v,\, \omega$")
    end

    % --- u
    subplot(313); hold on
    stairs(t_u,[u,u(end)],'LineWidth',linewidth,'LineStyle',stylek,'DisplayName',leg{k})
    if k == 1
        xlabel("$t$"); ylabel("$u$")
        grid on
        ylim([-1.1*u_max 1.1*u_max])
        yline(-u_max,'--','Color',matlab_blood_red,'HandleVisibility','off')
        yline( u_max,'--','Color',matlab_blood_red,'HandleVisibility','off')
    end
end

subplot(311); legend('Location','southeast');
subplot(312); legend('Location','southeast');
subplot(313); legend('Location','best');
