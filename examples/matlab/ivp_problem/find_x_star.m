syms x real
ts = -x/3;
L =(8/3*ts.^3 +8/3*x.*ts.^2 + 8/9*x.^2.*ts +1/3*T^3+1/3*x.*T^2+x.^2*T/9) + (T+(x-5)/3).^2;
nabla_L = diff(L,x)
x_star = solve(nabla_L)
a = double(x_star(1))
b = double(x_star(2))