clear all
import casadi.*
%%
settings = NosnocOptions();  
model = NosnocModel();  

settings.time_optimal_problem = 1;
settings.n_s = 2; 
settings.mpcc_mode = 'elastic_ineq';
settings.s_sot_max =	25;
settings.s_sot_min =	1/25;

settings.N_stages = 10; 
settings.N_finite_elements = 3; 
model.T = 1;    
q = SX.sym('q'); v = SX.sym('v'); 
model.x = [q;v]; model.x0 = [0;0]; 
model.lbx = [-inf;-25]; model.ubx = [inf;25];
u = SX.sym('u'); model.u = u;
model.lbu = -5; model.ubu = 5;
f_1 = [v;u]; f_2 = [v;3*u]; 
model.c = v-10; model.S = [-1;1]; model.F = [f_1 f_2];
model.g_terminal = [q-200;v-0];
solver = NosnocSolver(model, settings);
[results,stats] = solver.solve();

%% Plot
v_max = 25;
u_max = 5;
v_trash_hold = 10;
plot_results_nosnoc_tutorial
