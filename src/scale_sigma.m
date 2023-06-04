function out = scale_sigma(exprs, sigma, scale)
    import casadi.*
    symbolic_mode = ['casadi.' exprs.type_name()];
    out = define_casadi_symbolic(symbolic_mode, 'dummy', 0);
    for ii=1:size(exprs, 1)
        x = exprs(ii);
        out = vertcat(out, substitute(x, sigma, sigma/scale(ii)));
    end
end

