classdef MpecoptWrapper < handle & matlab.mixin.indexing.RedefinesParen
    properties
        mpcc
        opts
        stats
    end

    methods
        function obj = MpecoptWrapper(mpcc, opts)
        % Construct the wrapper for mpecopt. For now this doesn't do anything really :)
            obj.mpcc = mpcc;
            obj.opts = opts;
            obj.stats = struct;
        end

        function out = cat(dim,varargin)
            nosnoc.error('invalid', 'Invalid Operation')
        end

        function varargout = size(obj,varargin)
        % This is a required overload for matlab.mixin.indexing.RedefinesParen.
        % In the case of a scalar like this class returning 1 or throwing an error is prefered.
            varargout = {1};
        end
    end

    methods (Access=protected)
        function varargout = parenReference(obj, index_op)
            % TODO(@anton) this returns mpccresults struct from homotopy
            import casadi.*;
            p = inputParser;
            addParameter(p, 'x0', []);
            addParameter(p, 'lbx', []);
            addParameter(p, 'ubx', []);
            addParameter(p, 'lbg', []);
            addParameter(p, 'ubg', []);
            addParameter(p, 'p', []);
            addParameter(p, 'lam_g0', []);
            addParameter(p, 'lam_x0', []);
            parse(p, index_op(1).Indices{:});

            mpcc = obj.mpcc;
            opts = obj.opts;

            if isa(mpcc, 'vdx.problems.Mpcc')
                mpcc.finalize_assignments();

                if ~isempty(p.Results.x0)
                    mpcc.w.init = p.Results.x0;
                end
                if ~isempty(p.Results.lbx)
                    mpcc.w.lb = p.Results.lbx;
                end
                if ~isempty(p.Results.ubx)
                    mpcc.w.ub = p.Results.ubx;
                end
                if ~isempty(p.Results.lbg)
                    mpcc.g.lb = p.Results.lbg;
                end
                if ~isempty(p.Results.ubg)
                    mpcc.g.lb = p.Results.ubg;
                end
                if ~isempty(p.Results.p)
                    mpcc.p.val = p.Results.p;
                end
                if ~isempty(p.Results.lam_g0)
                    mpcc.g.init_mult = p.Results.lam_g0;
                end
                if ~isempty(p.Results.lam_x0)
                    mpcc.w.init_mult = p.Results.lam_x0;
                end

                mpcc_struct = mpcc.to_casadi_struct();
                solver_initialization = mpcc.to_solver_initialization();
                [solution, stats] = mpec_optimizer(mpcc_struct, solver_initialization, opts.mpecopt);
            else
                mpcc = obj.mpcc;
                opts = obj.opts;

                solver_intialization = p.Results;
                [solution, stats] = mpec_optimizer(mpcc, solver_initialization, opts.mpecopt)
            end
            varargout{1} = solution;
            obj.stats = stats;
        end

        function obj = parenAssign(obj,index_op,varargin)
            nosnoc.error('invalid', 'Invalid operation');
        end
        
        function obj = parenDelete(obj,index_op)
            nosnoc.error('invalid', 'Invalid operation')
        end

        function n = parenListLength(obj,index_op,ctx)
            n = 1;
        end
    end
end
