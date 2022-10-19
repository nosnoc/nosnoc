clear all
import casadi.*
%%
[settings] = default_settings_nosnoc();  
settings.time_optimal_problem = 1;
settings.n_s = 2; 
settings.mpcc_mode = 5;
settings.s_sot_max =	25;
settings.s_sot_min =	1/25;

model.N_stg = 10; model.N_FE = 3; model.T = 1;    
q = SX.sym('q'); v = SX.sym('v'); 
model.x = [q;v]; model.x0 = [0;0]; 
model.lbx = [-inf;-25]; model.ubx = [inf;25];
u = SX.sym('u'); model.u = u;
model.lbu = -5; model.ubu = 5;
f_1 = [v;u]; f_2 = [v;3*u]; 
model.c = v-10; model.S = [-1;1]; model.F = [f_1 f_2];
model.g_terminal = [q-200;v-0];
[results,stats,model,settings,solver_initialization] = nosnoc_solver(model,settings);

%% Plot
v_max = 25;
u_max = 5;
v_trash_hold = 10;
plot_results_nosnoc_tutorial
