function [model] = oscilator(model_in)

import casadi.*


%% Time horizon
if ~isempty(model_in)
unfold_struct(model_in,'caller');
else
    disp('empty struct input');
end
%% Model parameters

omega = 2*pi;
A1 = [1 omega;...
    -omega 1];
A2 = [1 -omega;...
    omega 1];

% T = 1+2*T_res;                            % time budget of transformed pseudo time
% T = 0.5+T_res;
% N_stages = 200;                   % number of control intervals
% N_stages = 2;                   % number of control intervals
% h = T/N_stages;                   % nominal step size (equidistant)


%% Inital Value
x0 = [exp(-1);0];
u0 = 0; % guess for control variables


%% Define model dimensions, equations, constraint functions, regions an so on.

n_simplex = 1;% number of Carteisna products in the model ("independet switches"), we call this layer
% number of modes in every simplex
m_1 = 2;
m_vec = [m_1];
%% Variable defintion
x1 = MX.sym('x1');
x2 = MX.sym('x2');
x = [x1;x2];
n_x = 2;

lbx = -inf*ones(n_x,1);
ubx = inf*ones(n_x,1);
% every constraint funcion corresponds to a simplex (note that the c_i might be vector valued)
c = x1^2+x2^2-1;
% sign matrix for the modes
S = [1;-1];
% discrimnant functions
h_1 = -S*c;

c = [c];
h_indictaros = [h_1];

%% modes of the ODEs layers (for all  i = 1,...,n_simplex);

% for c1, h1,
f_11 = A1*x;
f_12 = A2*x;


% in matrix form
f_1 = [f_11 f_12];
F = f_1;
% f = cat(3,f_1,f_2);


%% dummy control
u = MX.sym('u');
u0 = 0;
lbu = -1;
ubu = 1;
n_u = 1;
%% objective
f_q = 0*u^2+0;
f_q_T = 0;


%% Generic part
% (make of local workspace a struct and pass to output
names = who;

for ii = 1:length(names)
    eval([ 'model.' names{ii} '=' names{ii} ';'])
end

end

