% Solve an example nlp via FESD applied to the vector flow.
import casadi.*
x = SX.sym('x', 2);
f = (x(1)+1)^2 + (x(2))^2;
g = x(2)-x(1)^2 -2;
x0 = [1;5];
results = nlp_solver(x,f,g,x0)
