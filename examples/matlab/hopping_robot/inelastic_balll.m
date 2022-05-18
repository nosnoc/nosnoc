clear all
clear all;
import casadi.*
%%
[settings] = default_settings_nosnoc();  
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 2;
% settings.pss_mode = 'Step';
settings.pss_lift_step_functions= 0;
% settings.mpcc_mode = 2;
settings.N_homotopy = 3;
%%

model.N_stages = 10; % number of control intervals
model.N_finite_elements = 6; % number of finite element on every control intevral (optionally a vector might be passed)
model.T = 15;    % Time horizon
% Symbolic variables and bounds
a_n = 20;
q = SX.sym('q',2); 
v = SX.sym('v',2); 
t = SX.sym('t',1); 
model.x = [q;v;t]; % add all important data to the struct model,
model.x0 = [0;2;2;-1;0]; 
% control
% Dyanmics and the regions
f = [v;[0;-9.81];1];
f_aux_n = [0;0;0;a_n;0];
model.c = [q(2);v(2)];
if 1
    model.S = [1 1;...
              1 -1;...
           -1 1;...
           -1 -1];
    model.F = [f, f,f, f_aux_n];
else
model.S = [1 0;...
           -1 1;...
           -1 -1];
model.F = [f, f, f_aux_n];
end
%% Simulation setings
N_finite_elements = 3;
T_sim = 2;
N_sim = 20;

model.T_sim = T_sim;
model.N_finite_elements = N_finite_elements;
model.N_sim = N_sim;
settings.use_previous_solution_as_initial_guess = 1;
%% Call FESD Integrator
[results,stats,model] = integrator_fesd(model,settings);

%%
plot(results.x_res(1,:),results.x_res(2,:))