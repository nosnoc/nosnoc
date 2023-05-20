function [output] = solve_car_turbo_with_bonmin(model,settings)
%% Mixed integer problem with opti and bonmin.
import casadi.*
hybrid_dynamics = 1;
solve_relaxed = 0;
relax_terminal_constraint = 0;
time_optimal_problem = 1;
g_u_constraint = 0;
N_trails = settings.N_trails;
minlp_time_limit = 600;
%% collocation schme
% Degree of interpolating polynomial
n_s = settings.n_s;
% Get collocation points
collocation_scheme = settings.irk_scheme; % 'legendre' .  'radau'
tau = collocation_points(n_s, collocation_scheme);
% Collocation linear maps
[C,D,B] = collocation_coeff(tau);
%% paramseters
u_max = 5;
v_trash_hold = 10;
v_max = 25;
x_goal = [200 ;0];
%% Discretization
N_finite_elements = model.dims.N_finite_elements;
N_control_intevrlas = 1;
N_stages = model.dims.N_stages;

% Time horizon
T_val = 15; % Time horizon
opti = Opti();
T = opti.parameter();
h = T/(N_stages*N_control_intevrlas*N_finite_elements);
%

n_x = 2;
n_u = 1;
x0 = [0;0]; % Initial condition for x
u0 = [0];
% Declare model variables
q = MX.sym('q',1);
v = MX.sym('v',1);
u = MX.sym('u',n_u);
x = [q;v];

% Bounds on x
M = 1e3;
lbx = [-inf;-v_max];
ubx = [inf;v_max];
% Bounds on u
lbu = -u_max*ones(n_u,1);
ubu = u_max*ones(n_u,1);
% Model equations
f_A = [v;u];
f_B = [v;3*u];

% switching functions
psi = [v-v_trash_hold];
% Objective term
% if time_optimal_problem
L = 0;
% else
%     L = u'*u;
% end
% terminal costs
f_q_T = 0;
%% Formulate MINLP problem
bigM_minlp_formulation

%% solve with bonmin
cpu_time_all = [];
results = [];
options.time_limit = minlp_time_limit;
options.print_level = 1;
error_all = [];
sucess = 1;
discrete = struct('discrete', discrete);
opti.solver('bonmin',discrete,options);
% opti.solver('bonmin',discrete);
for kk = 1:N_trails
    tic
    try
        sol = opti.solve();
    catch
    end
    CPU_time_bonmin = opti.stats.t_wall_total;
    try
        x_opt = sol.value(Xs);
        u_opt = sol.value(Us);
        y_opt = sol.value(Ys);
        T_opt = sol.value(T_final);
    catch
        sucess = 0;
        T_opt = nan;
        %         x_opt = sol.value(Xs);
        %         u_opt = sol.value(Us);
        %         y_opt = sol.value(Ys);
        %         T_opt = sol.value(T_final);
    end
    CPU_time_bonmin1  = toc;
    if ~sucess
        CPU_time_bonmin = CPU_time_bonmin1;
    end
    fprintf('Tic toc bonmin time %2.2f s.\n',CPU_time_bonmin)
    cpu_time_all = [cpu_time_all,CPU_time_bonmin];
    if sucess
        [tout,yout,error] = car_turbo_sim(u_opt,T_opt,N_stages,hybrid_dynamics);
    else
        tout = nan;
        yout = nan;
        error = nan;
    end
    error_all = [error_all ;error];
end

if sucess
    results.x_opt = x_opt;
    results.u_opt = u_opt;
    results.T_opt = T_opt;
    results.y_opt = y_opt;
end


%%
fprintf('Total bonmin CPU time: %2.3f s.\n',CPU_time_bonmin);
if sucess
    [tout,yout,error] = car_turbo_sim(u_opt,T_opt,N_stages,hybrid_dynamics);
else
    tout = nan;
    yout = nan;
    error = nan;
end
%%
output.T_opt = T_opt;
output.error = error;
output.error_all = error_all;
output.tout = tout;
output.yout = yout;
output.results = results;
output.cpu_time_all = cpu_time_all;
output.cpu_time =  mean(cpu_time_all);
output.N_stages =  model.dims.N_stages;
output.N_finite_elements =  model.dims.N_finite_elements;

end


