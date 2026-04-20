% BSD 2-Clause License

% Copyright (c) 2022, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% This file is part of NOSNOC.

classdef Solver < handle & matlab.mixin.indexing.RedefinesParen
    properties
        mpcc % Either a struct with the (possibly optional) fields (f, p, w, g, G, H) or a subclass of vdx.problems.Mpcc.
        nlp
        opts % Options object.
        stats % Struct with stats of the last solve.
        solver
        casadi_opts
    end

    properties (Access=private)
        ind_scalar_G
        ind_scalar_H
        ind_nonscalar_G
        ind_nonscalar_H
        ind_map_G % Map from lifted indices to original indices in G.
        ind_map_H % Map from lifted indices to original indices in H.

        ind_cc1
        ind_cc2
        cc_types

        n_G_lift
        n_H_lift

        G_fun % Function (nlp.w, nlp.p)|-> mpcc.G
        H_fun % Function (nlp.w, nlp.p)|-> mpcc.G
        comp_res_fun % Function (nlp.w, nlp.p)|-> mmax(mpcc.G,mpcc.H)
        f_mpcc_fun % Function (nlp.w, nlp.p)|-> mpcc.f
        w_mpcc_fun % Function (nlp.w)|-> mpcc.w
        g_mpcc_fun % Function (nlp.w, nlp.p)|-> mpcc.g

        ind_map_g % Map for general constraints containing the corresponding indices in the original MPCC passed in and indices in the NLP
        ind_map_w % Map for primal variables containing the corresponding indices in the original MPCC passed in and indices in the NLP
        ind_map_p % Map for parameters containing the corresponding indices in the original MPCC passed in and indices in the NLP
    end

    methods(Access=public)

        function obj=Solver(mpcc, opts)
            import casadi.*
            if isa(mpcc, 'vdx.problems.Mpcc')
                mpcc = mpcc.to_casadi_struct();
            end
            casadi_symbolic_mode = split(class(mpcc.x), '.');
            casadi_symbolic_mode = casadi_symbolic_mode{end};
            % preprocess empty fields
            if ~isfield(mpcc, 'g') || isempty(mpcc.g)
                mpcc.g = casadi.(casadi_symbolic_mode)([]);
            end
            if ~isfield(mpcc, 'p') || isempty(mpcc.p)
                mpcc.p = casadi.(casadi_symbolic_mode)([]);
            end

            nx = length(mpcc.x);
            ng = length(mpcc.g);
            ncc = length(mpcc.G);

            ind_map_G = [];
            ind_map_H = [];
            [ind_scalar_G,ind_nonscalar_G, ind_map_G] = find_nonscalar(mpcc.G,mpcc.x,mpcc.p);
            [ind_scalar_H,ind_nonscalar_H, ind_map_H] = find_nonscalar(mpcc.H,mpcc.x,mpcc.p);
            obj.ind_nonscalar_G = ind_nonscalar_G;
            obj.ind_nonscalar_H = ind_nonscalar_H;
            obj.ind_scalar_G = ind_scalar_G;
            obj.ind_scalar_H = ind_scalar_H;
            obj.ind_map_G = ind_map_G;
            obj.ind_map_H = ind_map_H;

            n_G_lift = length(obj.ind_nonscalar_G);
            obj.n_G_lift = n_G_lift;
            n_H_lift = length(obj.ind_nonscalar_H);
            obj.n_H_lift = n_H_lift;

            % Based on what we know madmpec does. TODO(@anton) probably move this logic to casadi
            ind_cc1 = zeros(ncc,1);
            ind_cc2 = zeros(ncc,1);
            ind_cc1(ind_scalar_G) = ind_map_G;
            ind_cc2(ind_scalar_H) = ind_map_H;
            ind_cc1(ind_nonscalar_G) = (ng+1):(ng+n_G_lift);
            ind_cc2(ind_nonscalar_H) = (ng+1+n_G_lift):(ng+n_G_lift+n_H_lift);
            obj.ind_cc1 = ind_cc1;
            obj.ind_cc2 = ind_cc2;
            cctypes = zeros(ncc,1);
            cctypes(ind_nonscalar_G) = cctypes(ind_nonscalar_G) + 2;
            cctypes(ind_nonscalar_H) = cctypes(ind_nonscalar_H) + 1;

            casadi_opts = struct();
            casadi_opts.print_time = 0; % TODO(@anton) fix this
            casadi_opts.madnlp = opts.opts_madnlp;
            casadi_opts.ccopt = opts.opts_ccopt;
            casadi_opts.ind_cc = num2cell(horzcat(obj.ind_cc1, obj.ind_cc2),2);
            casadi_opts.cctypes = num2cell(cctypes);

            nlp = struct();
            nlp.f = mpcc.f;
            nlp.x = mpcc.x;
            nlp.g = vertcat(mpcc.g,mpcc.G(ind_nonscalar_G),mpcc.H(ind_nonscalar_H));
            nlp.p = mpcc.p;

            obj.mpcc = mpcc;
            obj.nlp = nlp;
            obj.opts = opts;

            obj.casadi_opts = casadi_opts;
            obj.solver = casadi.nlpsol('ccopt', 'ccopt', nlp, casadi_opts);

            G_fun = Function('G', {mpcc.x, mpcc.p}, {mpcc.G});
            H_fun = Function('H', {mpcc.x, mpcc.p}, {mpcc.H});
            obj.G_fun = G_fun;
            obj.H_fun = H_fun;
            comp_res_fun = Function('comp_res', {mpcc.x, mpcc.p}, {mmax(mpcc.G.*mpcc.H)});
            obj.comp_res_fun = comp_res_fun;
            f_mpcc_fun = Function('f_mpcc', {mpcc.x, mpcc.p}, {mpcc.f});
            obj.f_mpcc_fun = f_mpcc_fun;
        end

        function out = cat(dim,varargin)
            nosnoc.error('invalid', 'Invalid Operation')
        end

        function varargout = size(obj,varargin)
        % This is a required overload for matlab.mixin.indexing.RedefinesParen.
        % In the case of a scalar like this class returning 1 or throwing an error is prefered.
            varargout = {1};
        end

        function ind = end(obj,k,n)
        % This is a required overload for matlab.mixin.indexing.RedefinesParen.
        % In the case of a scalar like this class returning 1 or throwing an error is prefered.
            ind = 1;
        end
    end

    methods (Access=protected)
        function varargout = parenReference(obj, index_op)
            % TODO(@anton) this returns mpccresults struct from homotopy
            import casadi.*;
            p = inputParser;
            addParameter(p, 'x0', []);      % Initial guess.
            addParameter(p, 'y0', []);      % Initial complementarity active set guess (ignored by this solver).
            addParameter(p, 'lbx', []);     % Variable lower bounds.
            addParameter(p, 'ubx', []);     % Variable upper bounds.
            addParameter(p, 'lbg', []);     % Constraint lower bounds.
            addParameter(p, 'ubg', []);     % Constraint upper bounds.
            addParameter(p, 'p', []);       % Parameter upper bounds.
            addParameter(p, 'lam_g0', []);  % Initial constraint multipliers. (CasADi+ipopt convention)
            addParameter(p, 'lam_x0', []);  % Initial variable bound multipliers. (CasADi+ipopt convention)
            % TODO(@anton) Also perhaps initial complementarity multipliers but need to think about
            %              how exactly to implement that.
            parse(p, index_op(1).Indices{:});

            % data and local variables
            mpcc = obj.mpcc;
            nlp = obj.nlp;
            opts = obj.opts;

            nx = length(mpcc.x);
            nx_v = length(mpcc.x) + obj.n_G_lift + obj.n_H_lift;

            ng = length(mpcc.g);
            ng_v = length(mpcc.g) + length(obj.ind_nonscalar_G) + length(obj.ind_nonscalar_H);
            % Update nlp data
            np = length(mpcc.p);
            w_init = zeros(nx, 1);
            if ~isempty(p.Results.x0)
                w_init(1:nx) = p.Results.x0;
            end
            w_lb = zeros(nx, 1);
            if ~isempty(p.Results.lbx)
                w_lb(1:nx) = p.Results.lbx;
            end
            w_ub = zeros(nx, 1);
            if ~isempty(p.Results.ubx)
                w_ub(1:nx) = p.Results.ubx;
            end
            g_lb = zeros(ng_v, 1);
            if ~isempty(p.Results.lbg)
                g_lb(1:ng) = p.Results.lbg;
            end
            g_ub = zeros(ng_v, 1);
            g_ub(ng+1:end) = inf;
            if ~isempty(p.Results.ubg)
                g_ub(1:ng) = p.Results.ubg;
            end
            p_val = zeros(np, 1);
            if ~isempty(p.Results.p)
                p_val = p.Results.p;
            end
            g_init_mult = zeros(ng_v,1);
            if ~isempty(p.Results.lam_g0)
                g_init_mult(1:ng) = p.Results.lam_g0;
            end
            w_init_mult = zeros(nx, 1);
            if ~isempty(p.Results.lam_x0)
                w_init_mult(1:nx) = p.Results.lam_x0;
            end

            if any(w_lb(obj.ind_map_G) == -inf)
                ind_add_lb_G = find(w_lb(obj.ind_map_G) == -inf);
                warning("adding lb=0 to x1 vars that are missing them")
                w_lb(obj.ind_map_G(ind_add_lb_G)) = 0.0;
            end

            if any(w_lb(obj.ind_map_H) == -inf)
                ind_add_lb_H = find(w_lb(obj.ind_map_H) == -inf);
                warning("adding lb=0 to x2 vars that are missing them")
                w_lb(obj.ind_map_H(ind_add_lb_H)) = 0.0;
            end

            nlp_results = obj.solver('x0', w_init,...
                'lbx', w_lb,...
                'ubx', w_ub,...
                'lbg', g_lb,...
                'ubg', g_ub,...
                'lam_g0', g_init_mult,...% TODO(@anton) perhaps we use init instead of mult.
                'lam_x0', w_init_mult,...
                'p', p_val);

            mpcc_results = nlp_results;
            mpcc_results.f = full(nlp_results.f);
            mpcc_results.x = full(nlp_results.x(1:nx));
            mpcc_results.g = full(nlp_results.g(1:ng));
            mpcc_results.G = full(obj.G_fun(mpcc_results.x, p_val));
            mpcc_results.H = full(obj.H_fun(mpcc_results.x, p_val));
            mpcc_results.lam_g = full(nlp_results.lam_g(1:ng));
            mpcc_results.lam_x = full(nlp_results.lam_x(1:nx)); % TODO(@anton) probably set x1 x2 lam_x to zero
            mpcc_results.lam_G = full(obj.solver.stats.ccopt.multipliers_x1);
            mpcc_results.lam_H = full(obj.solver.stats.ccopt.multipliers_x2);

            varargout{1} = mpcc_results;
            obj.stats = obj.solver.stats;
            obj.stats.converged = obj.stats.success;
            obj.stats.constraint_violation = max(obj.stats.ccopt.primal_feas,...
                obj.stats.ccopt.cc_feas);
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

        function print_infeasibility(obj, results)
            if obj.opts.print_details_if_infeasible
                print_problem_details(results,obj.model,obj.mpcc, []);
            end
            if obj.opts.pause_homotopy_solver_if_infeasible
                keyboard
            end
        end
    end
end


function X = all_combinations(varargin)
    numSets = length(varargin);
    for i=1:numSets,
        thisSet = sort(varargin{i});
        if ~isequal(prod(size(thisSet)),length(thisSet)),
            nosnoc.error('combinations_not_vectors', 'All inputs must be vectors.')
        end
        if ~isnumeric(thisSet),
            nosnoc.error('cominations_not_numeric','All inputs must be numeric.')
        end
        sizeThisSet(i) = length(thisSet);
        varargin{i} = thisSet;
    end
    X = zeros(prod(sizeThisSet),numSets);
    for i=1:size(X,1)
        ixVect = cell(length(sizeThisSet),1);
        sz = flip(sizeThisSet);
        if length(sz) == 1
            sz = [sz,1];
        end
        [ixVect{:}] = ind2sub(sz,i);
        ixVect = flip([ixVect{:}]);
        vect = zeros(1, numSets);
        for jj=1:numSets
            vect(jj) = varargin{jj}(ixVect(jj));
        end
        X(i,:) = vect;
    end
end
