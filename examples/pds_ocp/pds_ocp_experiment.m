clear all
close all
import casadi.*

%% Options
n_s = 3;
T = 10;
x_target = [0;-4];

Ns = [5, 10, 15, 20, 25, 50, 100];
%% Define and solve OCP
fesd_errors = [];
h_ks = [];
for N_stages=Ns
    [x_res, u_res, h_k] = solve_pds_ocp(true, N_stages, n_s, T, x_target);
    x_sim = integrate_system(N_stages, T, u_res);
    fesd_errors = [fesd_errors, norm(x_sim(:,end)-x_target)];
    h_ks = [h_ks, h_k];
end
no_errors = [];
for N_stages=Ns
    [x_res, u_res, h_k] = solve_pds_ocp(false, N_stages, n_s, T, x_target);
    x_sim = integrate_system(N_stages, T, u_res);
    no_errors = [no_errors, norm(x_sim(:,end)-x_target)];
end
%% plot results
fontsize = 18;

figure
hold on
plot(h_ks,no_errors(1)*h_ks.^2,':k', "LineWidth", 2, "DisplayName", "$O(h^2)$")
plot(h_ks,fesd_errors(1)*h_ks.^5,'--k', "LineWidth", 2, "DisplayName", "$O(h^5)$")
plot(h_ks, fesd_errors, '-x', "LineWidth", 3, "Color", [27,158,119]/256, "MarkerSize", 14, "DisplayName", "FESD");
plot(h_ks, no_errors, '-x', "LineWidth", 3, "Color", [217,95,2]/256, "MarkerSize", 14, "DisplayName", "Standard");
ylabel("$x(T)-\bar{x}$")
xlabel("$\bar{h}$")
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
ax = gca;
ax.FontSize = fontsize;

legend("Location", "southeast")
hold off
