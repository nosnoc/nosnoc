classdef Pss < nosnoc.model.Base
% A model that represents a Piecewise Smooth System with regions $R_i\subset \R^{n_x}$ and dynamics $f_i(u,x)$.
% $$\xdot = f_i(x,u),\ \mathrm{if}\ x\in R_i$$
% for $i= 1,\ldots,n_f$.
%
% This system is convexified via the Filippov convexification which under an existance assumptions $\dot{x}$ is:
% $$\xdot\in\overline{\conv}\Set{f_i(x,u)}{i = 1,\ldots,n_f}$$
%
% Todo:
%     Document independent systems.
%
% See Also:
%    More information about this model and its discretization can be found in :cite:`Nurkanovic2024a`.
    properties
        F % cell<casadi.SX|casadi.MX>: $n_x \times n_f$ matrix whose columns define the vector fields in the correspoding PSS regions $R_i$.
        S % cell<casadi.SX|casadi.MX>: Sign matrix.
        c % cell<casadi.SX|casadi.MX>: Switching functions defining region boundaries.
        g_ind % cell<casadi.SX|casadi.MX>: Vector of length $n_f$ containing all Stewart indicator functions
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

            if isempty(obj.g_ind)
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
                    obj.g_ind{ii} = -obj.S{ii}*obj.c{ii};
                    % dimensions of c
                    dims.n_c_sys  = [dims.n_c_sys;length(obj.c{ii})];
                end
            else
                if ~iscell(obj.g_ind)
                    obj.g_ind = {obj.g_ind};
                end
            end

            dims.n_f_sys = arrayfun(@(sys) size(obj.F{sys},2),1:dims.n_sys);

            obj.dims = dims;
        end
    end
end
