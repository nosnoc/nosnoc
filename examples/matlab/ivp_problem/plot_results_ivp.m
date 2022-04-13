
%% read solutions

 unfold_struct(results,'base');

%% plots
figure
plot(t_grid,x_opt);
hold on
xlabel('$t$','interpreter','latex');
ylabel('$x(t)$','interpreter','latex');
grid on
figure
semilogy(stats.complementarity_stats,'k')
xlabel('iter','interpreter','latex');
ylabel('complementarity','interpreter','latex');
grid on
