classdef Pds < nosnoc.model.Base
    properties
        f_x
        c
        E
    end

    methods
        function obj = Pds()
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

            if size(obj.f_x,1) ~= dims.n_x
                error("nosnoc: f_x has incorrect dimension. It must have the same dimension as x.")
            end

            dims.n_c = size(obj.c,1);

            if ~isempty(E)
                if size(E, 1) ~= size(E,2)
                    error("nosnoc: Projection matrix E must be square.");
                end
                if size(E, 1) ~= dims.n_c
                    error("nosnoc: Projection matrix E must be an n_c by n_c matrix where n_c is the number of functions defining the set C.")
                end
            else
                E = eye(dims.n_c);
            end

            obj.dims = dims;
        end
    end
end
