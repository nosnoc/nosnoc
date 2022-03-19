function [U,V] = oscilator_eval(X,Y)
%OSCILATOR_EVAL Summary of this function goes here
%   Detailed explanation goes here
omega = 2*pi;
A1 = [1 omega;...
    -omega 1];
A2 = [1 -omega;...
    omega 1];
U = [];
V = [];
for ii=1:size(X,1)
    x = X(ii,:);
    y = Y(ii,:);
    u = 0.5*(1+sign(x.^2+y.^2-1)).*(A1(1,1).*x+A1(1,2).*y) + 0.5*(1-sign(x.^2+y.^2-1)).*(A2(1,1).*x+A2(1,2).*y);
    v = 0.5*(1+sign(x.^2+y.^2-1)).*(A1(2,1).*x+A1(2,2).*y) + 0.5*(1-sign(x.^2+y.^2-1)).*(A2(2,1).*x+A2(2,2).*y);
    U = [U;u];
    V = [V;v];
end


end

