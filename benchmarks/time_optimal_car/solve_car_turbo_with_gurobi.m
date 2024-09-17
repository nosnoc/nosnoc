function [output] = solve_car_turbo_with_gurobi(model,settings)
%% Mixed integer problem with opti gurobi
import casadi.*
hybrid_dynamics = 1;
solve_relaxed = 0;
relax_terminal_constraint = 0;
time_optimal_problem = 0;
solve_gurobi = 1;
use_Q_matrix = 0;
g_u_constraint = 0;
N_trails = settings.N_trails;
verbose_gurobi  = 1;
plot_gurobi = 0;
%% collocation schme
% Degree of interpolating polynomial
n_s = settings.n_s;
% Get collocation points
collocation_scheme = settings.rk_scheme; % 'legendre' .  'radau'
tau = collocation_points(n_s, collocation_scheme);
% Collocation linear maps
[C,D,B] = collocation_coeff(tau);

%% Discretization
N_finite_elements = settings.N_finite_elements;
N_control_intevrlas = 1;
N_stages = settings.N_stages;

% Time horizon
T_val = 15; % Time horizon
opti = Opti();
T = opti.parameter();
h = T/(N_stages*N_control_intevrlas*N_finite_elements);
%%
u_max = 5;
v_trash_hold = 10;
v_max = 25;
x_goal = [200 ;0];
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
M = 1e6;
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
if time_optimal_problem
    L = 0;
else
    L = u'*u;
end
% terminal costs
f_q_T = 0;
%% Formulate MINLP problem
bigM_minlp_formulation
error_all = [];
cpu_time_all = [];
%% model functions
for kk = 1:N_trails
    CPU_time_gurobi = nan;
    T_ub_k = 20;
    T_lb_k = 1;
    N_iter = 20;
    iter = nan;
    tol_biseciton = 1e-6;

    CPU_time_gurobi = nan;
    if solve_gurobi
        gurobi_bisection_solver
        T_final_opt_gurobi = T_final_opt;
    end
     cpu_time_all = [cpu_time_all,CPU_time_gurobi];
    % Plot gurobi results
    % Read and plot results
    %     w_opt_gurobi = result.x;
    u_opt = w_opt_gurobi(ind_u);
    u_opt = reshape(u_opt,n_u,length(u_opt)/n_u);
    x_opt = w_opt_gurobi(ind_x);
    x_opt = reshape(x_opt,n_x,length(x_opt)/n_x);
    y_opt = w_opt_gurobi(ind_y);
    y_opt = reshape(y_opt,n_y,length(y_opt)/n_y);
    q_opt = x_opt(1,:);
    v_opt = x_opt(2,:);
    T_opt = T_val_k(end);
    [tout,yout,error] = car_turbo_sim(u_opt,T_final_opt,N_stages,hybrid_dynamics);
    error_all = [error_all;error];

end
results.x_opt = x_opt;
results.u_opt = u_opt;
results.T_opt = T_opt;
results.y_opt = y_opt;

%%
fprintf('Total gurobi CPU time: %2.3f s.\n',CPU_time_gurobi);
[tout,yout,error] = car_turbo_sim(u_opt,T_final_opt,N_stages,hybrid_dynamics);
%%
output.T_opt = T_final_opt;
output.error = error;
output.tout = tout;
output.yout = yout;
output.error = error;
output.error_all = error_all;
output.results = results;
output.cpu_time_all = cpu_time_all;
output.cpu_time =  mean(cpu_time_all);
output.N_stages =  settings.N_stages;
output.N_finite_elements =  settings.N_finite_elements ;
output.iter=  iter;






end


