clear all
clc
close all
import casadi.*
%% Settings
[settings] = default_settings_fesd();  %% Optionally call this function to have an overview of all options.
settings.time_freezing = 1; settings.d = 3; 
model.T = 5; model.N_stages = 15; model.N_finite_elements = 3;
model.x0 = [0;0.5;0;0;0];
q = MX.sym('q',2); v = MX.sym('v',2); t = MX.sym('t');
model.x = [q;v;t];
model.c = q(2); model.S = [1; -1];
u = MX.sym('u',2); model.u = u;
f_1 = [v;u(1);u(2)-9.81;1]; f_2 = [0;v(2);0;-10*(q(2))-0.211989*v(2);0];
model.F = [f_1 f_2];
model.f_q = u'*u; model.f_q_T = 10*v'*v;
model.g_ineq = u'*u-7^2;
model.g_terminal = q-[4;0.25];
[solver,solver_initalization, model,settings] = create_nlp_fesd(model,settings);
[results,stats] = homotopy_solver(solver,model,settings,solver_initalization);
plot_results_thowring_ball
