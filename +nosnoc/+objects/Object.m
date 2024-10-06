classdef Object < handle & matlab.mixin.Heterogeneous
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
