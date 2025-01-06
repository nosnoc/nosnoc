classdef Counter < handle
% A class designed for keeping a class-static counter similar to the `static` keyword in a c++ class.
    properties
        count = 0;
    end

    properties(Constant)
        counter = nosnoc.objects.Counter();
    end
    
    methods(Access=private)
        function obj = Counter()

        end
    end
end 
