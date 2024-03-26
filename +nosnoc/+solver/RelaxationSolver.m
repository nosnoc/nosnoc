classdef RelaxationSolver < handle & matlab.mixin.indexing.RedefinesParen
    properties
        mpcc
        nlp
        opts
        stats % struct
    end

    methods (Access=public)
        
        function obj=RelaxationSolver(mpcc, opts)
        % TODO(@anton) move NosnocSolver building here.
            import casadi.*
            casadi_symbolic_mode = class(mpcc.w);
            obj.mpcc = mpcc;
            obj.opts = opts;
            opts.preprocess();
            stats = struct;
            
            if isa(mpcc, 'vdx.problems.Mpcc') % We can properly interleave complementarities if we get a vdx.Mpcc
                use_vdx = true;
                % TODO(@anton) implement this
            else % Otherwise use vdx internally anyway but be sad about interleaving
                use_vdx = false;
                nlp = vdx.Problem('casadi_type', casadi_symbolic_mode);
                nlp.w.mpcc_w(0) = {mpcc.w};
                nlp.p.mpcc_p(0) = {mpcc.p};
                nlp.g.mpcc_g(0) = {mpcc.g};
                nlp.p.sigma_p(0) = {{'sigma_p', 1}, 0, inf, opts.sigma_0};

                switch opts.elasticity_mode
                  case ElasticityMode.NONE
                    sigma = nlp.p.sigma_p(0);
                  case ElasticityMode.ELL_INF
                    nlp.w.s_elastic(0) = {{'s_elastic', 1}, opts.s_elastic_min, opts.s_elastic_max, s_elastic_0};
                    sigma = nlp.w.s_elastic(0);
                  case ElasticityMode.ELL_1
                    n_c = size(mpcc.G);
                    nlp.w.s_elastic(0) = {{'s_elastic', n_c}, opts.s_elastic_min, opts.s_elastic_max, s_elastic_0};
                    sigma = nlp.w.s_elastic(0);
                end

                % add complementarites
                expr = psi_fun(mpcc.G, mpcc.H, sigma);
                [lb, ub, expr] = generate_mpcc_relaxation_bounds(expr, opts);
                nlp.g.complementarities(0) = {expr, lb, ub};

                if ~opts.assume_lower_bounds % Lower bounds on G, H, not already present in MPCC
                   nlp.g.G_lower_bounds(0) = {mpcc.G, 0, inf};
                   nlp.g.H_lower_bounds(0) = {mpcc.G, 0, inf};
                end
            end
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
            addParameter(p, 'lbx', 'SX');
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
