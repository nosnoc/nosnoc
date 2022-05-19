clear all;
close all;
clc;
import casadi.*
%% settings
[settings] = default_settings_nosnoc();  
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 2;

settings.pss_mode = 'Step';
settings.pss_lift_step_functions= 1;

settings.print_level = 4;
settings.mpcc_mode = 3;
settings.opts_ipopt.ipopt.max_iter = 5e2;
settings.N_homotopy = 12;
settings.sigma_0 = 10;
settings.cross_comp_mode = 10;

settings.initial_lambda_0 = 0.5; settings.initial_lambda_1 = 0.5; settings.initial_alpha = 0.5;

settings.use_fesd = 1;
settings.time_freezing = 1;


settings.step_equilibration = 1;
settings.step_equilibration_mode = 3;
%%
% Symbolic variables and bounds
q = SX.sym('q',2); v = SX.sym('v',2); 
u = SX.sym('u',1);  u_max = 10;
model.x = [q;v]; 
model.u = u;
model.lbu = -u_max;
model.ubu = u_max;
model.e = 0;
model.mu = 0.25;
model.a_n = 20;
model.f = [u;-9.81];
%%
a1 = 1;
a2 = 1.2;
a3 = -1;
sigma_c = 1e-1;
model.c = q(2)-(a1*sin(a2*q(1))+a3*q(1))*(1-tanh((q(1))/sigma_c))/2;

qx0 = -15;
qy0 = (a1*sin(a2*qx0)+a3*qx0 )*(1-tanh((qx0 )/sigma_c))/2;
qx_target = 4;
qy_target = (a1*sin(a2*qx_target )+a3*qx_target )*(1-tanh((qx_target)/sigma_c))/2;
model.x0 = [qx0;qy0;0;0]; 
model.f_q = u^2;
model.g_terminal = q-[qx_target;qy_target];
%% Discretization
model.T= 5;
model.N_stages = 50;
model.N_FE = 3;
%% NOSNOC OCP solver
[results,stats,model] = nosnoc_solver(model,settings);
%%
qx = results.x_opt(1,:);
qy = results.x_opt(2,:);
vx = results.x_opt(3,:);
vy = results.x_opt(4,:);
t_opt = results.x_opt(5,:);
%%
h  = model.h/2;
figure
tt = -20:0.05:25;
for ii = 1:length(qx);
plot(tt,(a1*sin(a2*tt)+a3*tt).*(1-tanh((tt)/sigma_c))/2,'Color',0.2*ones(3,1),'LineWidth',1)
hold on
axis equal
plot(qx(ii),qy(ii),'LineWidth',2,'Marker','.','MarkerSize',15);
grid on
pause(h)
hold off
clf
end

%%
figure
plot(t_opt,vy);
hold on
plot(t_opt,vx);
%%
