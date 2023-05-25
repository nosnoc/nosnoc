function plot_two_ball_traj(results, title_str)

if nargin < 2
    title_str = '';
end
benchmark_globals;
N_FE = NFE_VALUES(1);

%% read and plot results
q1 = results.x(1,:);
q2 = results.x(2,:);
v1 = results.x(3,:);
v2 = results.x(4,:);
t_grid = results.t_grid;

%%
figure
subplot(311)
plot(t_grid,q1,'LineWidth',1.5);
hold on
plot(t_grid,q2,'LineWidth',1.5);
yline(R,'k--')
xlim([0 t_grid(end)])
% ylim([-1.0 max([q1,q2])+1])
grid on
ylabel('$q$','interpreter','latex');
xlabel('$t$','interpreter','latex');
% axis equal
subplot(312)
plot(t_grid,v1,'LineWidth',1.5);
hold on
plot(t_grid,v2,'LineWidth',1.5);
xlim([0 t_grid(end)])
ylim([-max(abs([v1,v2]))-1.0 max(abs([v1,v2]))+1])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
subplot(313)
Lambda_opt = [results.Lambda_normal];
stem(t_grid(1:N_FE:end),[Lambda_opt,nan])
hold on
xlim([-0.01 t_grid(end)])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda_{\mathrm{n}}$','interpreter','latex');
title(title_str)

end