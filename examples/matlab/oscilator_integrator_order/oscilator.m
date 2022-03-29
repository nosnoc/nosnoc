function [model] = oscilator(model_in)

import casadi.*
%% Time horizon
if ~isempty(model_in)
    unfold_struct(model_in,'caller');
else
    disp('empty struct input');
    smooth_model = 0;
end
%% Model parameters

omega = 2*pi;
A1 = [1 omega;...
    -omega 1];
A2 = [1 -omega;...
    omega 1];
if smooth_model
    A2 = A1;
end
%% Inital Value
x0 = [exp(-1);0];
%% Variable defintion
x1 = MX.sym('x1');
x2 = MX.sym('x2');
x = [x1;x2];
% every constraint funcion corresponds to a simplex (note that the c_i might be vector valued)
c = x1^2+x2^2-1;
% sign matrix for the modes
S = [1;-1];
c = [c];

f_11 = A1*x;
f_12 = A2*x;
% in matrix form
F = [f_11 f_12];
%% Generic part
% (make of local workspace a struct and pass to output
names = who;

for ii = 1:length(names)
    eval([ 'model.' names{ii} '=' names{ii} ';'])
end

end

