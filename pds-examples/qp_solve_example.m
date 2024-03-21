% Solve an example QP via FESD applied to the vector flow.
import casadi.*
if 0
    Q = -2*eye(2);
    c = [2;2];
    A = [eye(2);-1,-1];
    b = [0;0;-4];
    x0 = [1.1;1.05];
else
    n = 20;
    m = 10;
    Q = rand(n,n);
    Q = 10*Q'*Q;
    c = 0*rand(n,1);
    A = rand(m,n);
    b = rand(m,1);

    % start at a feasbile point
    if 0 % optimization approach
        x = SX.sym('x', n);

        fqp = struct('x', x, 'f', 0, 'g', A*x - b);
        fqp_sol = nlpsol('fqp', 'ipopt', fqp);
        fqp_res = fqp_sol('lbg',0);

        x0 = full(fqp_res.x);
    else % just do some linear algebra approach
        x0 = A\(b+0.1);
    end
end
[results,prob] = qp_solver(Q,c,A,b,x0)
