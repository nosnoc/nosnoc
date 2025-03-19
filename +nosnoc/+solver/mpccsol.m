function mpccsol = mpccsol(name, plugin, mpcc, options)
% MPCCSOL create an mpccsolver functor in the style of casadi nlpsol
% TODO(@anton) If we want to make this an actual casADi function we probably need to implement the solvers in C++.
%              For now we just return a class with an overriden paren

    if strcmp(plugin, 'reg_homotopy')
        mpccsol = nosnoc.reg_homotopy.Solver(mpcc, options);
    elseif strcmp(plugin, 'mpecopt')
        mpccsol = mpecopt.Solver(mpcc, options);
    end
end
