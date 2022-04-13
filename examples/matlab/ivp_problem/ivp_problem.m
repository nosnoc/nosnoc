function [model] = ivp_problem(varargin)

import casadi.*
if nargin > 0
    unfold_struct(varargin{1},'caller')
else
    %% Time horizon
T = 2;
N_stages = 1;
N_finite_elements = 25;
x0 = -1;
end

%% Variable defintion
x = MX.sym('x');
c = x;
S = [1;-1];
%% modes of the ODEs layers (for all  i = 1,...,n_simplex);
f_11 = [1];
f_12 = [3];
F = [f_11 f_12];
%% objective
f_q = x^2;
f_q_T = (x-5/3)^2;
%% Generic part

names = who;
for ii = 1:length(names)
    eval([ 'model.' names{ii} '=' names{ii} ';'])
end

end

