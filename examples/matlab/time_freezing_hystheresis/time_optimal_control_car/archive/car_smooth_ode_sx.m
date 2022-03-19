close all;
clear all
import casadi.*


%  Colors
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
organe = [0.9290 0.6940 0.1250];
grey = [0.85 0.85 0.85];


%%
T = 12;
q_goal = 150;
v_goal = 0;
u_max = 5;

time_optimal_problem = 1;
sot_max = 10;

%% Inital Value

q0 = 0;
v0 = 0;
L0 = 0;
a0 = 0;
t0 = 0;


x0 = [q0;v0;L0;a0;t0];

%%
% Degree of interpolating polynomial
d = 3;
% Get collocation points
tau = collocation_points(d, 'radau');
% Collocation linear maps
[C,D,B] = collocation_coeff(tau);
% Time horizon
% Control discretization
N = 100; % number of control intervals
h = T/N;
% Declare model variables
q = SX.sym('q');
v = SX.sym('v');
L = SX.sym('L');
a = SX.sym('a');
t = SX.sym('t');
x = [q;v;L;a;t];
u1 = SX.sym('u1');
sot = SX.sym('sot');
if time_optimal_problem
    u = [u1;sot];
    else
    u = [u1];
end

n_x = 5;
n_u = length(u);

Pn = 10;
% Model equations
if time_optimal_problem
    xdot = sot*[v;u1;Pn;1;1];
else
    xdot = [v;u;Pn;1;1];
end


if time_optimal_problem
    obj = sot;
else
    obj = u1^2;
end

% Continuous time dynamics
f = Function('f', {x, u}, {xdot, obj});

% Start with an empty NLP

opti = Opti();
J = 0;

u_max = 5;


if time_optimal_problem
    lbu = [-u_max;sot_max^(-1)];
    ubu = [u_max; sot_max];
    u0 = [0 T/N];
else
    lbu = -u_max;
    ubu = u_max;
    u0 = 0;
end


% "Lift" initial conditions
Xk = opti.variable(n_x);
opti.subject_to(Xk==x0);
opti.set_initial(Xk, x0);

% Collect all states/controls
Xs = {Xk};
Us = {};

% Formulate the NLP
for k=0:N-1
    % New NLP variable for the control
    Uk = opti.variable(n_u);
    Us{end+1} = Uk;
    opti.subject_to(lbu <= Uk);
    opti.subject_to(Uk <=ubu);
    opti.set_initial(Uk, u0);
    
    % Decision variables for helper states at each collocation point
    Xc = opti.variable(n_x, d);
    %    opti.subject_to(-0.25 <= Xc(1,:));
    opti.set_initial(Xc, repmat(x0,1,d));
    
    % Evaluate ODE right-hand-side at all helper states
    [ode, quad] = f(Xc, Uk);
    
    % Add contribution to quadrature function
    J = J + quad*B*h;
    
    % Get interpolating points of collocation polynomial
    Z = [Xk Xc];
    
    % Get slope of interpolating polynomial (normalized)
    Pidot = Z*C;
    % Match with ODE right-hand-side
    opti.subject_to(Pidot == h*ode);
    
    % State at end of collocation interval
    Xk_end = Z*D;
    
    % New decision variable for state at end of interval
    Xk = opti.variable(n_x);
    Xs{end+1} = Xk;
    opti.set_initial(Xk, x0);
    
    % Continuity constraints
    opti.subject_to(Xk_end==Xk)
end
opti.subject_to(Xk_end(1:2)==[q_goal;v_goal])

Xs = [Xs{:}];
Us = [Us{:}];

opti.minimize(J);

opti.solver('ipopt');

cpu_times_all = [];
for ii = 1:5
tic
sol = opti.solve();
cpu_time = toc;
cpu_times_all = [cpu_times_all;cpu_time];
end


x_opt = sol.value(Xs);
u_opt = sol.value(Us);

u_opt1 = u_opt(1,:);
if time_optimal_problem
    sot_opt = u_opt(2,:);
end

%% Plot the solution
tgrid = linspace(0, T, N+1);
q_opt = x_opt(1,:);
v_opt = x_opt(2,:);
L_opt = x_opt(3,:);
a_opt = x_opt(4,:);
t_opt = x_opt(5,:);

%%
figure
tgrid = 0:h:T;
subplot(221)
plot(t_opt,q_opt);
xlabel('$t$','Interpreter','latex');
ylabel('$q(t)$','Interpreter','latex');
subplot(222)
plot(t_opt,v_opt);
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$v(t)$','Interpreter','latex');
subplot(223)
plot(t_opt,t_opt);
hold on
grid on
if time_optimal_problem
stairs(t_opt, [sot_opt(1:end) nan])
end
xlabel('$t$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
subplot(224)
stairs(t_opt, [u_opt1(1:end) nan])
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$u(t)$','Interpreter','latex');


T_final = t_opt(end)
v_max = max(v_opt)
