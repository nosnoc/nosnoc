classdef pss < nosnoc.model.base
    properties
        F
        S
        c
        g_ind
    end

    methods
        function obj = pss()
            obj = obj@nosnoc.model.base();
        end
        
        function verify_and_backfill(obj, opts)
            arguments
                obj
                opts nosnoc.Options
            end
            import casadi.*
            verify_and_backfill@nosnoc.model.base(obj,opts);

            dims = obj.dims;

            % check how many subsystems are present
            if iscell(obj.F)
                dims.n_sys = length(obj.F);
            else
                obj.F = {obj.F};
                dims.n_sys = 1;
            end
            
            % Check if all data is avilable and if dimensions match.
            if ~iscell(obj.S)
                obj.S = {obj.S};
            end
            if length(obj.S) ~= dims.n_sys
                error('nosnoc: Number of matrices S does not match number of subsystems. Note that the number of subsystems is taken to be number of matrices F_i which collect the modes of every subsystem.')
            end
            % Check constraint function c
            if isempty(obj.c)
                error('nosnoc: Expresion for c, the constraint function for regions R_i is not provided.');
            else
                if ~iscell(obj.c)
                    obj.c = {obj.c};
                end
                if length(obj.c) ~= dims.n_sys
                    error('nosnoc: Number of different expressions for c does not match number of subsystems (taken to be number of matrices F_i which collect the modes of every subsystem).')
                end
            end

            dims.n_c_sys = [];
            for ii = 1:dims.n_sys
                if size(obj.S{ii},2) ~= length(obj.c{ii})
                    error('nosnoc: The matrix S and vector c do not have compatible dimension.');
                end
                
                % dimensions of c
                dims.n_c_sys  = [dims.n_c_sys;length(obj.c{ii})];
            end

            dims.n_f_sys = arrayfun(@(sys) size(obj.F{sys},2),1:dims.n_sys);

            obj.dims = dims;
        end
    end
end
