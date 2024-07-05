clear all
close all
import casadi.*

%% Options
use_fesd = true;
N_stages = 4;
n_s = 3;
T = 10;
x_target = [0;-4];

%% Define and solve OCP
[x_res, u_res] = solve_pds_ocp(use_fesd, N_stages, n_s, T, x_target);

%% plot OCP results
fontsize = 12;

figure
hold on
% unperturbed vector field
[X,Y] = meshgrid(-6:1:6, -6:1:6);
U = -0.2*(X+1).^2;
V = -0.4*(Y+3);
quiver(X,Y,U,V,'r')
plot(x_res(1,:), x_res(2,:), '-b', "LineWidth", 3)

xlim([-6,6])
ylim([-6,6])

a=2;
b=4; 
x0=0;
y0=0;
t=-pi:0.01:pi;
x=x0+a*cos(t);
y=y0+b*sin(t);
plot(x,y, '--r', "LineWidth", 2)
axis equal

%% Make integrator and plot results
x_sim = integrate_system(N_stages, T, u_res);
plot(x_sim(1,:), x_sim(2,:), '-g', "LineWidth", 3)

error = norm(x_sim(:,end)-x_target);

disp(['Terminal error e=', num2str(error), ' for n_s=', num2str(n_s), ' N_stages=', num2str(N_stages), ' and use_fesd=', mat2str(use_fesd)])
