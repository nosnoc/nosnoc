function mpccsol = mpccsol(name, plugin, mpcc, options)
% MPCCSOL create an mpccsolver functor in the style of casadi nlpsol
% TODO(@anton) If we want to make this an actual casADi function we probably need to implement the solvers in C++.
%              For now we just return a class with an overriden paren

    mpccsol = nosnoc.solver.RelaxationSolver(plugin, mpcc, options);
end
