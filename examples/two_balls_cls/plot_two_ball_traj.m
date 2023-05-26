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
figure('Renderer', 'painters', 'Position', [100 100 900 220])
subplot(121)
plot(t_grid,q1,'LineWidth',1.5);
hold on
plot(t_grid,q2,'LineWidth',1.5);
yline(R,'k--')
xlim([0 t_grid(end)])
ylim([-0.1 2])
grid on
yline(0.2,'k-','LineWidth',1.5)
ylabel('$q$','interpreter','latex');
xlabel('$t$ ','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
legend('ball 1', 'ball 2', 'interpreter','latex')
% axis equal
subplot(122)
plot(t_grid,v1,'LineWidth',1.5);
hold on
plot(t_grid,v2,'LineWidth',1.5);
xlim([0 t_grid(end)])
ylim([-max(abs([v1,v2]))-1.0 max(abs([v1,v2]))+1])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
legend('ball 1', 'ball 2', 'interpreter','latex')

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
eval(['print -dpdf -painters ' ['two_balls_states'] ])

%%
figure
Lambda_opt = [results.Lambda_normal];
stem(t_grid(1:N_FE:end),[Lambda_opt,nan],'LineWidth',1.5);
hold on
xlim([-0.01 t_grid(end)])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda_{\mathrm{n}}$','interpreter','latex');
title(title_str)

end