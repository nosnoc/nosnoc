classdef base < matlab.mixin.Scalar
    properties
    end

    methods(Abstract)
        generate_variables(obj, opts)
        generate_equations(obj, opts)
    end
end
