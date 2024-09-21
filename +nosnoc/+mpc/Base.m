classdef (Abstract) Base < matlab.mixin.Scalar & handle
    methods (Abstract)
        [u, stats] = get_feedback(obj, x0);
        [stats] = do_preparation(obj);
    end
end

