%% Mixed integer problem with opti and bonmin.
clear all
import casadi.*

hybrid_dynamics = 1;
solve_relaxed = 0;
relax_terminal_constraint = 0;
time_optimal_problem = 1;
solve_bonmin = 0;
solve_gurobi = 1;
use_Q_matrix = 0;

if solve_gurobi
    time_optimal_problem = 0; % no nonlinearites, make grid search
end
%%
verbose_gurobi = 1;
plot_gurobi = 1;
%% collocation schme
% Degree of interpolating polynomial
n_s = 2;
% Get collocation points
collocation_scheme = 'radau'; % 'legendre' .  'radau'
tau = collocation_points(n_s, collocation_scheme);
% Collocation linear maps
[C,D,B] = collocation_coeff(tau);

%% Discretization
N_finite_elements = 3;
N_control_intevrlas = 1;
N_stages = 10;

% Time horizon
T_val  = 15;
opti = Opti();
T = opti.parameter();
h = T/(N_stages*N_control_intevrlas*N_finite_elements);
X_goal = [150;0];
%% Problem defintion
n_x = 2;
n_u = 1;

% Declare model variables
q = MX.sym('q');
v = MX.sym('v');
x = [q;v];
u = MX.sym('u',n_u);
u_max = 5;

% Model parameters
x0 = [0;0]; % Initial condition for x
u0 = 0;

v_max = 25;
v_1 = 15;
v_2 = 10;
% Bounds on x
lbx = [0;-v_max];
ubx = [inf;v_max];

% Bounds on u
lbu = -u_max*ones(n_u,1);
ubu = u_max*ones(n_u,1);
% Model equations
f_A = [v;u];
f_B = [v;3*u];
% Transition functions
l_AB = [v_1-v];
l_BA = [v-v_2];

% l_AB = [v_2-v];
% l_BA = [v-v_1];

T_final_opt_bonmin = nan;
T_final_opt_gurobi = nan; 

% Objective term
if time_optimal_problem
    L = 0;
else
    L = 1*u^2;
end
%% Formulate MINLP
formulate_costas_minlp
%% Solve via Bonmin/Ipopt
if solve_relaxed
    opti.solver('ipopt');
else
    opti.solver('bonmin',struct('discrete', discrete));
end

CPU_time_bonmin  = 0;
if solve_bonmin
    sol = opti.solve();
    CPU_time_bonmin = opti.stats.t_wall_total;

    x_opt = sol.value(Xs);
    u_opt = sol.value(Us);
    y_opt = sol.value(Ys);
    L_AB_n_opt = sol.value(L_transition_n_s);
    L_AB_opt = sol.value(L_transition_s);
    if time_optimal_problem
        T_final_opt = sol.value(T_final);
    else
        T_final_opt = 1;
    end
    T_final_opt_bonmin = T_final_opt;
    fprintf('Objective value : %2.3f \n',full(sol.value(J)));
    % Plot the solution
    plot_results_hystersis_car      
end
%% gurobi bisection solve
T_ub_k = 20;
T_lb_k = 5;
N_iter = 20;
iter = nan;
tol_biseciton = 1e-12;

CPU_time_gurobi = nan;
if solve_gurobi
   gurobi_bisection_solver
   plot_results_hystersis_car  
   T_final_opt_gurobi = T_final_opt;
end
% w_opt_gurobi_lp = w_opt_gurobi;

%% empirical
% u_opt_long = repmat(u_opt,N_finite_elements,1);
% u_opt_long= u_opt_long(:)';
% h = T_final_opt/(N_finite_elements*N_stages*N_control_intevrlas);
% emprical_u = (diff(v_opt(1:end))/(h))./u_opt_long;
% %%
%%
fprintf('Final time gurobi: %2.6f \n',T_final_opt_gurobi);
fprintf('Final time bonmin: %2.6f \n',T_final_opt_bonmin);
fprintf('CPU time bonmin: %2.3f s.\n',CPU_time_bonmin);
fprintf(' CPU time gurobi : %2.3f s.\n',CPU_time_gurobi);
fprintf('Bisection iters: %d.\n',iter);
%% check lp/qp
 
[tout,yout,error] = car_hysteresis_sim(u_opt,T_final_opt,N_stages);

% use_Q_matrix = 0;
% T_val = 12.14;
% use_Q_matrix = 1;
% miqp_set_up_problem
% model.obj = ones(n_w,1);
% % rmfield(model,'Q')
% model.start = w_opt_gurobi_lp;
% result = gurobi(model, params)
% % %
% w_opt_gurobi = result.x;
% u_opt = w_opt_gurobi(ind_u);
% u_opt = reshape(u_opt,n_u,length(u_opt)/n_u);
% x_opt = w_opt_gurobi(ind_x);
% x_opt = reshape(x_opt,n_x,length(x_opt)/n_x);
% y_opt = w_opt_gurobi(ind_y);
% y_opt = reshape(y_opt,n_y,length(y_opt)/n_y);
% plot_results_hystersis_car