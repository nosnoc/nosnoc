clear all
clc
close all
import casadi.*
%% Settings
[settings] = default_settings_fesd();  %% Optionally call this function to have an overview of all options.
% FESD

%% Time settings
% settings.d = 3;
settings.time_freezing = 1;
%% Generate Model
% model = throwing_ball_model();
model.T = 5; model.N_stages = 15; model.N_finite_elements = 3;
model.q0 = [0;0.5];
qf = [6;0.25];
model.x0 = [0;0.5;0;0;0];
% model.h = T/N_stages;
k_tf = 10;  gamma_tf = 0.9; % restitution coefficient
c_tf = 2*abs(log(gamma_tf))*((k_tf) /(pi^2+log(gamma_tf)^2) )^(1/2);
T_res = 2*pi/sqrt((4*k_tf-c_tf^2));
g = 9.811;

%% Define model dimensions, equations, constraint functions, regions an so on.
n_simplex = 1;
m_1 = 2;
m_vec = [m_1];
%% Variable defintion
qx = MX.sym('qx'); qy = MX.sym('qy');
vx = MX.sym('vx'); vy = MX.sym('vy');
t = MX.sym('t');
q = [qx;qy];v = [vx;vy]; x = [q;v;t];
n_x = length(x);
n_q = 2;
lbx = -inf*ones(n_x,1);
ubx = inf*ones(n_x,1);
% every constraint funcion corresponds to a simplex (note that the c_i might be vector valued)
b_bottom = 0;
c_1 = qy-b_bottom; % bottom
S1 = [1; -1];
% discrimnant functions
h_1 = -S1*[c_1];
c = [c_1];
h_indictaros = [h_1];
% control
ux = MX.sym('ux'); uy = MX.sym('uy'); u = [ux;uy];
n_u = 2;
u0 = [0;0];
umax = inf;
if active_control
    lbu = -umax*ones(n_u,1);
    ubu = umax*ones(n_u,1);
else
    lbu = 0*ones(n_u,1);
    ubu = 0*ones(n_u,1);
end

%% modes of the ODEs layers (for all  i = 1,...,n_simplex);
% for c1, h1,

% f_aux_top = -[0;vy;0;-k_tf*qy-c_tf*vy;0];
% f_aux_left = [vx;0;-k_tf*qx-c_tf*vx;0;0];

f_aux_bottom = [0;vy;0;-k_tf*(qy-b_bottom)-c_tf*vy;0];

f_11 = [vx;vy;ux;uy-g;1];
f_12 = f_aux_bottom;
% in matrix form
f_1 = [f_11 f_12];

%% objective


f_q = u'*u;
% Terminal Cost
f_q_T = 1*v'*v;

%%  general nonlinear constinrst
general_nonlinear_constraint  = 1;
% g_ineq = v'*v;
g_ineq = u'*u;
g_ineq_lb = [-inf];
g_ineq_ub = [v_max_R^2];
g_ineq_fun  = Function('g_ineq_fun',{x,u},{g_ineq});
% g_ineq_fun  = Function('g_ineq_fun',{x},{g_ineq});
%% Terminal constraint
terminal_constraint = 1;
g_terminal = q-qf;
g_terminal_lb = zeros(2,1);
g_terminal_ub = zeros(2,1);
g_terminal_fun = Function('g_terminal_fun',{x},{g_terminal});
%% Reformulation of the PSS into a DCS
[model,settings] = model_reformulation_fesd(model,settings);
%% Formulate NLP;
if 1
    [solver,solver_initalization, model,settings] = create_nlp_fesd(model,settings);
else
    [solver,solver_initalization, model,settings] = create_nlp_mfe_develop(model,settings);
end
%% Get variables into main workspace
unfold_struct(model,'base');
unfold_struct(settings,'base');
unfold_struct(solver_initalization,'base');
%% Solve NLP
[sol,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);
obj = full(sol.f);
w_opt = full(sol.x);
%% Read and plot Result
plot_results_thowring_ball
