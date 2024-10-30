classdef Object < handle & matlab.mixin.Heterogeneous
% Base class for nosnoc objects:
% - Ball (2 or 3 dimensions)
% - (union of) Ellipsoid (2 or 3 dimensions)
    properties(GetAccess=public, SetAccess=protected)
        n_dim
        name
    end

    properties(Constant)
        num_objects = nosnoc.objects.Counter();
    end

    methods (Sealed=true)
        function is_eq = eq(a,b)
            is_eq = eq@handle(a,b);
        end
    end
end
