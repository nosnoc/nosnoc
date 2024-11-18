classdef Pss < nosnoc.model.Base
    properties
        F % n_x times n_f matrix whose columns define the vector fields in the correspoding PSS regions
        S % sign matrix
        c % switching functions defining region boundaries
        g_indicator % Vector of length n_f containing all Stewart indicator functions
    end

    methods
        function obj = Pss()
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

            if isempty(obj.g_indicator)
                if length(obj.S) ~= dims.n_sys
                    nosnoc.error('size_S', 'Number of matrices S does not match number of subsystems. Note that the number of subsystems is taken to be number of matrices F_i which collect the modes of every subsystem.')
                end
                % Check constraint function c
                if isempty(obj.c)
                    nosnoc.error('missing_c', 'Expresion for c, the constraint function for regions R_i is not provided.');
                else
                    if ~iscell(obj.c)
                        obj.c = {obj.c};
                    end
                    if length(obj.c) ~= dims.n_sys
                        nosnoc.error('size_c', 'Number of different expressions for c does not match number of subsystems (taken to be number of matrices F_i which collect the modes of every subsystem).')
                    end
                end

                dims.n_c_sys = [];
                for ii = 1:dims.n_sys
                    if size(obj.S{ii},2) ~= length(obj.c{ii})
                        nosnoc.error('c_S_mismatch', 'The matrix S and vector c do not have compatible dimension.');
                    end
                    obj.g_indicator{ii} = -obj.S{ii}*obj.c{ii};
                    % dimensions of c
                    dims.n_c_sys  = [dims.n_c_sys;length(obj.c{ii})];
                end
            else
                if ~iscell(obj.g_indicator)
                    obj.g_indicator = {obj.g_indicator};
                end
            end

            dims.n_f_sys = arrayfun(@(sys) size(obj.F{sys},2),1:dims.n_sys);

            obj.dims = dims;
        end
    end
end
