classdef Object < handle & matlab.mixin.Heterogeneous
% Base class for nosnoc objects:
% - Ball (2 or 3 dimensions)
% - (union of) Ellipsoid (2 or 3 dimensions)
    properties(GetAccess=public, SetAccess=protected)
        n_dim % Number of dimensions (2 or 3).
        name % Name used in variable creation.
    end

    properties(Constant)
        num_objects = nosnoc.objects.Counter.counter; % Class static counter for number of objects of a given type.
    end

    properties(Access=public)
        timestamp
    end

    methods (Sealed=true)
        function is_eq = eq(a,b)
        % Check for equality via handle equality.
            is_eq = eq@handle(a,b);
        end
    end
end
