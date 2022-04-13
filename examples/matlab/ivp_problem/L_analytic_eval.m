function [L_analytic] = L_analytic_eval(x0)
T = 2;
ts = -x0/3;
L_analytic =(8/3*ts.^3 +8/3*x0.*ts.^2 + 8/9*x0.^2.*ts +1/3*T^3+1/3*x0.*T^2+x0.^2*T/9) + (T+(x0-5)/3).^2;

end