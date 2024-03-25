classdef RelaxationSolver < handle & matlab.mixin.indexing.RedefinesParen
    properties
        mpcc
        opts
        stats % struct
    end

    methods (Access=public)
        
        function obj=RelaxationSolver(mpcc, opts)
        % TODO(@anton) move NosnocSolver building here.
        end

        function out = cat(dim,varargin)
            error('Concatenation not supported')
        end

        function varargout = size(obj,varargin)
            varargout = 1;
            %TODO(anton) needs to return correct values for varargin
        end
        
        function ind = end(obj,k,n)
            ind = 1;
        end
    end

    methods (Access=protected)
        function varargout = parenReference(obj, index_op)
        % TODO(@anton) this returns mpccresults struct from homotopy
            p = inputParser;
            addParameter(p, 'casadi_type', 'SX');
            parse(p, varargin{:});
        end
        
        function obj = parenAssign(obj,index_op,varargin)
            error('Invalid operation');
        end
        
        function obj = parenDelete(obj,index_op)
            error('Invalid operation')
        end

        function n = parenListLength(obj,index_op,ctx)
            n = 1;
        end
    end
end
