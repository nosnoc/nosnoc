classdef Heaviside < nosnoc.model.Base
    properties
        f_x

        c
        alpha
    end

    methods
        function obj = Heaviside()
            obj = obj@nosnoc.model.Base();
        end
        
        function verify_and_backfill(obj, opts)
            arguments
                obj
                opts nosnoc.Options
            end
            import casadi.*
            verify_and_backfill@nosnoc.model.Base(obj,opts);

            dims = obj.dims;

            % check how many subsystems are present
            if iscell(obj.c)
                dims.n_sys = length(obj.c);
            else
                obj.c = {obj.c};
                dims.n_sys = 1;
            end

            dims.n_c_sys = [];
            for ii = 1:dims.n_sys
                % dimensions of c
                dims.n_c_sys  = [dims.n_c_sys;length(obj.c{ii})];
            end

            if size(obj.alpha,1) ~= sum(dims.n_c_sys)
                error('nosnoc: There needs to be a step function alpha for each switching function in c');
            end
            dims.n_alpha = size(obj.alpha, 1);

            obj.dims = dims;
        end
    end
end
