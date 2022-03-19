function [model] = temp_control_model_voronoi()
import casadi.*
%% Discretization parameters
N_stages = 2;
N_finite_elements = 1;
T = 0.1; % (here determined latter depeding on omega)
h = T/N_stages;
%% Model Parameters
active_control = 0;
% inital value
t0 = 0;
w0 = 0;
y0 = 15;

lambda_cool_down = -0.2; % cool down time constant of lin dynamics
u_heat = 10; % heater power

% jump points in x in the hysteresis function
y1 = 18;
y2 = 20;


z1 = [1/4;-1/4];
z2 = [1/4;1/4];
z3 = [3/4;3/4];
z4 = [3/4;5/4];
% Z = [1/4 1/4 3/4 3/4;...
%      -1/4 1/4 3/4 5/4]
Z = [z1 z2 z3 z4];

% if 0
% % Z = [x1 x2 x3 x4; w1 w2 w3 w4];
% x = linspace(-2,2,100);
% i = 1; j = 2;
% y12 = -((Z(1,i)-Z(1,j))./(Z(2,i)-Z(2,j)))*x+0.5*(norm(Z(:,i))^2-norm(Z(:,j))^2)/((Z(2,i)-Z(2,j)));
% i = 2; j = 3;
% y23 = -((Z(1,i)-Z(1,j))./(Z(2,i)-Z(2,j)))*x+0.5*(norm(Z(:,i))^2-norm(Z(:,j))^2)/((Z(2,i)-Z(2,j)));
% i = 3; j = 4;
% y34 = -((Z(1,i)-Z(1,j))./(Z(2,i)-Z(2,j)))*x+0.5*(norm(Z(:,i))^2-norm(Z(:,j))^2)/((Z(2,i)-Z(2,j)));
% i = 2; j = 4;
% y24 = -((Z(1,i)-Z(1,j))./(Z(2,i)-Z(2,j)))*x+0.5*(norm(Z(:,i))^2-norm(Z(:,j))^2)/((Z(2,i)-Z(2,j)));
% i = 1; j = 3;
% y13 = -((Z(1,i)-Z(1,j))./(Z(2,i)-Z(2,j)))*x+0.5*(norm(Z(:,i))^2-norm(Z(:,j))^2)/((Z(2,i)-Z(2,j)));
% figure
% plot(x,y12)
% hold on
% grid on
% plot(x,y23)
% plot(x,y34)
% plot(Z(1,:),Z(2,:),'ko')
% end

%% Inital Value
x0 = [y0;w0;t0];

%% Model parameters for time freezing

%% Define model dimensions, equations, constraint functions, regions an so on.
n_simplex = 1;% number of Carteisna products in the model ("independet switches"), we call this layer
m_1 = 4;
m_vec = [m_1];
%% Variable defintion
y = MX.sym('y');
w = MX.sym('w');
t = MX.sym('t');

x = [y;w;t];
n_x = length(x);
lbx = -inf*ones(n_x,1);
ubx = inf*ones(n_x,1);
% every constraint funcion corresponds to a simplex (note that the c_i might be vector valued)
% c = [1;1;1];

% linear transformation for rescaling of the switching function.
psi = (y-y1)/(y2-y1);
z = [psi;w];
% discriminant functions via voronoi
h_1 = -2*z'*z1+norm(z1)^2;
h_2 = -2*z'*z2+norm(z2)^2;
h_3 = -2*z'*z3+norm(z3)^2;
h_4 = -2*z'*z4+norm(z4)^2;

% 
h_11 = norm([psi;w]-z1)^2;
h_12 = norm([psi;w]-z2)^2;
h_13 = norm([psi;w]-z3)^2;
h_14 = norm([psi;w]-z4)^2;

% h_1 = 1;
% h_2 = 1;
% h_3 = 15;
% h_4 = 25;

h_1 = [h_11;h_12;h_13;h_14];
h_indictaros = [h_1];
c = h_indictaros;
% h_indictaros = [1;1;15;15];

%% control
u = MX.sym('u');
n_u = 1;
u0 = [0];

if active_control
    umax = 1e-3;
    lbu = -umax*ones(n_u,1);
    ubu = umax*ones(n_u,1);
else
    lbu = 0*ones(n_u,1);
    ubu = 0*ones(n_u,1);
end

%% modes of the ODEs layers (for all  i = 1,...,n_simplex);
%
u_heat = 10;
f_A = [lambda_cool_down*y+u_heat;0;1];
f_B = [lambda_cool_down*y;0;1];

a_push = 5;
f_push_down = [0;-a_push*(psi-1)^2/(1+(psi-1)^2);0];
f_push_up = [0;a_push*(psi)^2/(1+(psi)^2);0];

f_11 = 2*f_A-f_push_down;
f_12 = f_push_down;
f_13 = f_push_up;
f_14 = 2*f_B-f_push_up;

% f_11 = ones(3,1);
% f_12 = ones(3,1);
% f_13 = ones(3,1);
% f_14 = ones(3,1);



f_1 = [f_11 f_12 f_13 f_14];

%% objective
f_q = active_control*(u^2)+y^2;
% Terminal Cost
f_q_T = 0;

%%  general nonlinear constinrst
general_nonlinear_constraint  = 0;
g_ineq = u^2;
g_ineq_lb = [-inf];
g_ineq_ub = [inf];
g_ineq_fun  = Function('g_ineq_fun',{x,u},{g_ineq});
%% Generic part
% (make of local workspace a struct and pass to output
names = who;
for ii = 1:length(names)
    eval([ 'model.' names{ii} '=' names{ii} ';'])
end
end

