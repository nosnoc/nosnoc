opti.set_value(T, T_val)
lbg = full(evalf(opti.lbg));
ubg = full(evalf(opti.ubg));

b0 = sparse(g_fun(w0,T_val));
A_dense = sparse(jac_g_fun(w0,T_val));

obj = full(nabla_f_fun(w0,T_val));
Q = sparse(hess_f_fun(w0,T_val));

A  = [];
b = [];
sense = [];
for ii = 1:length(lbg)
    if lbg(ii) == ubg(ii)
        A = [A;A_dense(ii,:)];
        sense = [sense,'='];
        b = [b;-b0(ii)+lbg(ii)];
    else
        if ~isequal(lbg(ii),-inf)
            A = [A;A_dense(ii,:)];
            sense = [sense,'>'];
            b = [b;-b0(ii)+lbg(ii)];
        end
        if ~isequal(ubg(ii),inf)
            A = [A;A_dense(ii,:)];
            sense = [sense,'<'];
            b = [b;-b0(ii)+ubg(ii)];
        end
    end
end
% solve with gurobi
model.A = A;

if use_Q_matrix
    model.Q  = Q;
else
    try
        rmfield(model,"Q");
    catch
    end

end
if ~use_Q_matrix
    if length(unique(obj)) == 1
        if (unique(obj)) == 00
            obj = ones(n_w,1);
        end
    end
end
model.obj = obj;
model.rhs = b;
model.sense = sense;
model.vtype = vtype;
model.lb = -inf*ones(n_w,1);
model.modelsense = 'min';
params.outputflag = 0;