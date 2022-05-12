import casadi.*
    n_w = length(discrete);
    w = opti.x;
    w0 = zeros(n_w,1);
    g = opti.g;
    f = opti.f;
    nabla_f = f.jacobian(w)';
    hess_f  = f.hessian(w);
    jac_g = g.jacobian(w);

    g_fun = Function('g_fun',{w,T},{g});
    jac_g_fun = Function('jac_g_fun',{w,T},{jac_g});
    f_fun = Function('f_fun',{w},{f});
    nabla_f_fun = Function('nabla_f_fun',{w,T},{nabla_f});
    hess_f_fun = Function('hess_f_fun',{w,T},{hess_f});

    % set up variable typs for miqp
       % MIQP set up
    
    vtype  = [];
    for ii = 1:length(discrete)
        if discrete(ii) == 1
            vtype  = [vtype 'B'];
        else
            vtype  = [vtype 'C'];
        end
    end