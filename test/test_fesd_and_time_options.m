function [results,stats,model,settings] = git
%TEST_SIMPLE_CAR_MODEL Test the simple car model accross cross
%complementarity and mpcc modes
import casadi.*
[settings] = default_settings_nosnoc();  
settings.irk_scheme = 'Radau-IIA';
% set the cross complimentarity mode that 
settings.cross_comp_mode = 3;
settings.mpcc_mode = 3;
settings.n_s = 2;
settings.homotopy_parameter_rule = 'superlinear';
settings.N_homotopy = 8;
% settings being tested
settings.use_fesd = use_fesd;
settings.time_optimal_problem = time_optimal_problem;
settings.equidistant_control_grid = equidistant_control_grid;
settings.use_speed_of_time_variables = use_speed_of_time_variables;
settings.local_speed_of_time_variable = local_speed_of_time_variable;

% Model - define all problem functions and
% Discretization parameters
model.N_stages = 10; % number of control intervals
model.N_finite_elements = 3; % number of finite element on every control intevral (optionally a vector might be passed)
% Symbolic variables and bounds
q = SX.sym('q'); v = SX.sym('v'); 
model.x = [q;v]; % add all important data to the struct model,
model.x0 = [0;0]; % inital value
% bounds on states
model.lbx = [-inf;-20];
model.ubx = [inf;20];
% control
u = SX.sym('u'); model.u = u;
model.lbu = -5; model.ubu = 5;
% Dyanmics and the regions
f_1 = [v;u]; % mode 1 - nominal
f_2 = [v;3*u]; % mode 2 - turbo
model.S = [-1;1];
model.F = [f_1 f_2];
model.c = v-10;
% Add terminal constraint
model.g_terminal = [q-200;v-0];
if time_optimal_problem
    model.T = 1;    
else
    model.T = 20;   
    model.f_q = u^2;
end
% Solve OCP
[results,stats,model,settings] = nosnoc_solver(model,settings);
if ~isempty(results.T_opt) && results.T_opt < 1e-2
    warning('Something went wrong.')
    disp(results.T_opt)
end
disp(results.f_opt)

end

