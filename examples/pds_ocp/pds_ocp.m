clear;
close all
import casadi.*

%% Options
use_fesd = true;
N_stages = 10;
n_s = 3;
T = 10;
x_target = [0;-4];

%% Define and solve OCP
[x_res, u_res, t_res, t_control, lambda_res, c_res] = solve_pds_ocp(use_fesd, N_stages, n_s, T, x_target);

%% plot OCP results
fontsize = 18;

figure
hold on
% unperturbed vector field
[X,Y] = meshgrid(-6:1:6, -6:1:6);
U = -0.2*(X+1).^2;
V = -0.4*(Y+3);
quiver(X,Y,U,V, "DisplayName", "$f(x)$", "Color", [27,158,119]/256, "LineWidth", 1.5)
a=2;
b=4; 
x0=0;
y0=0;
t=-pi:0.01:pi;
x=x0+a*cos(t);
y=y0+b*sin(t);
plot(x,y, '--r', "LineWidth", 2, "DisplayName", "$c(x)$")
plot(x_res(1,:), x_res(2,:), '-b', "LineWidth", 3, "DisplayName", "$x(t)$")
xlim([-5.5,5.5])
ylim([-5.5,5.5])
ax = gca;
ax.FontSize = fontsize;
axis equal
xlabel("$x$")
ylabel("$y$")
legend("Location", "northeastoutside")
hold off
grid on

% Plot controls
figure
stairs(t_control, [u_res, u_res(:, end)]', "LineWidth", 2)
xlabel("$t$")
ylabel("$u(t)$")
ax = gca;
ax.FontSize = fontsize;
ylim([-1.1, 1.1])
legend(["$u_x(t)$", "$u_y(t)$"])
grid on


figure
subplot(2,1,1)
plot(t_res, lambda_res, "LineWidth", 2)
xlabel("$t$")
ylabel("$\lambda(t)$")
ax = gca;
ax.FontSize = fontsize;
subplot(2,1,2)
plot(t_res, c_res, "LineWidth", 2)
xlabel("$t$")
ylabel("$c(x(t))$")
grid on
%% Make integrator and plot results
x_sim = integrate_system(N_stages, T, u_res);

error = norm(x_sim(:,end)-x_target);

disp(['Terminal error e=', num2str(error), ' for n_s=', num2str(n_s), ' N_stages=', num2str(N_stages), ' and use_fesd=', mat2str(use_fesd)])
