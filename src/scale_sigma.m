function out = scale_sigma(exprs, sigma, scale)
    import casadi.*
    out = [];
    for ii=1:size(exprs, 1)
        x = exprs(ii);
        out = vertcat(out, substitute(x, sigma, sigma/scale(ii)));
    end
end

