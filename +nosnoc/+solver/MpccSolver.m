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

% TODO(@anton) add printing

classdef MpccSolver < handle & matlab.mixin.indexing.RedefinesParen
    properties
        mpcc
        nlp
        opts
        stats % struct
        nlp_solver
        plugin
        relaxation_type
    end

    properties (Access=private)
        ind_scalar_G
        ind_scalar_H
        ind_nonscalar_G
        ind_nonscalar_H
        ind_map_G
        ind_map_H
    end

    methods (Access=public)
        
        function obj=MpccSolver(relaxation_type, mpcc, opts)
            import casadi.*
            import nosnoc.solver.*
            casadi_symbolic_mode = split(class(mpcc.x), '.');
            casadi_symbolic_mode = casadi_symbolic_mode{end};
            % preprocess empty fields
            if ~isfield(mpcc, 'g') || isempty(mpcc.g)
                mpcc.g = casadi.(casadi_symbolic_mode)([]); 
            end
            if ~isfield(mpcc, 'p') || isempty(mpcc.p)
                mpcc.p = casadi.(casadi_symbolic_mode)([]); 
            end
            obj.relaxation_type = relaxation_type;
            obj.mpcc = mpcc;
            obj.opts = opts;
            opts.preprocess();
            stats = struct;
            
            if isa(mpcc, 'vdx.problems.Mpcc') % We can properly interleave complementarities if we get a vdx.Mpcc
                use_vdx = true;
                % TODO(@anton) implement this
            else % Otherwise use vdx internally anyway but be sad about interleaving
                use_vdx = false;

                % create an nlp (in the form of a vdx.Problem) from the mpcc data (f,x,p,g).
                % We track the indices of this data in vdx.Variables: mpcc_w, mpcc_p, mpcc_g.
                nlp = vdx.Problem('casadi_type', casadi_symbolic_mode);
                nlp.w.mpcc_w = {mpcc.x}; % nlp.w.mpcc_w is a vdx.Variable (depth 0) which stores mpcc optimization variables
                nlp.p.mpcc_p = {mpcc.p}; % nlp.p.mpcc_p is a vdx.Variable (depth 0) which stores mpcc parameters
                nlp.g.mpcc_g = {mpcc.g}; % nlp.g.mpcc_g is a vdx.Variable (depth 0) which stores mpcc general constraints
                nlp.p.sigma_p = {{'sigma_p', 1}, opts.sigma_0}; % sigma_p is
                nlp.f = mpcc.f;
                n_c = length(mpcc.G);

                % Create relaxation slacks/parameters
                switch opts.homotopy_steering_strategy
                  case HomotopySteeringStrategy.DIRECT
                    % nlp.p.sigma_p(): sigma is a parameter/variable that has no indices
                    sigma = nlp.p.sigma_p(); 
                  case HomotopySteeringStrategy.ELL_INF
                    % adding elastic variables to nlp.w which augments the original mpcc.w
                    nlp.w.s_elastic = {{'s_elastic', 1}, opts.s_elastic_min, opts.s_elastic_max, opts.s_elastic_0};
                    sigma = nlp.w.s_elastic();
                    if opts.objective_scaling_direct
                        nlp.f = nlp.f + (1/nlp.p.sigma_p())*sigma;
                    else
                        nlp.f = nlp.p.sigma_p()*nlp.f + sigma;
                    end
                  case HomotopySteeringStrategy.ELL_1
                    % adding elastic variables to nlp.w which augments the original mpcc.w
                    nlp.w.s_elastic = {{'s_elastic', n_c}, opts.s_elastic_min, opts.s_elastic_max, opts.s_elastic_0};
                    sigma = nlp.w.s_elastic();
                    sum_elastic = sum1(sigma);
                    if opts.objective_scaling_direct
                        nlp.f = nlp.f + (1/nlp.p.sigma_p())*sum_elastic;
                    else
                        nlp.f = nlp.p.sigma_p()*nlp.f + sum_elastic;
                    end
                end

                [ind_scalar_G,ind_nonscalar_G, ind_map_G] = find_nonscalar(mpcc.G,mpcc.x);
                [ind_scalar_H,ind_nonscalar_H, ind_map_H] = find_nonscalar(mpcc.H,mpcc.x);
                obj.ind_nonscalar_G = ind_nonscalar_G;
                obj.ind_nonscalar_H = ind_nonscalar_H;
                obj.ind_scalar_G = ind_scalar_G;
                obj.ind_scalar_H = ind_scalar_H;
                obj.ind_map_G = ind_map_G;
                obj.ind_map_H = ind_map_H;
                % possibly lift complementarities
                if opts.lift_complementarities
                    ind_G_fun = Function('ind_G', {mpcc.x}, {mpcc.G.jacobian(mpcc.x)});
                    [ind_G1,ind_G2] = find(sparse(DM(ind_G_fun(nlp.w.sym) == 1)));
                    ind_H_fun = Function('ind_G', {mpcc.x}, {mpcc.H.jacobian(mpcc.x)});
                    [ind_H1,ind_H2] = find(sparse(DM(ind_H_fun(nlp.w.sym) == 1)));

                    nlp.w.G_lift = {{'G', length(ind_nonscalar_G)}, 0, inf};
                    G = casadi.(casadi_symbolic_mode)(size(mpcc.H,1), 1);
                    G(ind_scalar_G) = mpcc.G(ind_scalar_G);
                    G(ind_nonscalar_G) = nlp.w.G_lift();

                    
                    nlp.w.H_lift = {{'H', length(ind_nonscalar_H)}, 0, inf};
                    H = casadi.(casadi_symbolic_mode)(size(mpcc.H,1), 1);
                    H(ind_scalar_H) = mpcc.H(ind_scalar_H);
                    H(ind_nonscalar_H) = nlp.w.H_lift();
                    
                    nlp.g.G_lift = {mpcc.G(ind_nonscalar_G)-G(ind_nonscalar_G)};
                    nlp.g.H_lift = {mpcc.H(ind_nonscalar_H)-H(ind_nonscalar_H)};
                    
                else
                    G = mpcc.G;
                    H = mpcc.H;
                end

                % apply relaxation
                psi_fun = get_psi_fun(MpccMethod(obj.relaxation_type), opts.normalize_homotopy_update);
                lb = [];
                ub = [];
                expr = [];
                n_comp_pairs = size(G, 1);
                for ii=1:n_comp_pairs
                    expr_i = psi_fun(G(ii), H(ii), sigma);
                    [lb_i, ub_i, expr_i] = generate_mpcc_relaxation_bounds(expr_i, obj.relaxation_type);
                    lb = [lb;lb_i];
                    ub = [ub;ub_i];
                    expr = [expr;expr_i];
                end
                nlp.g.complementarities(0) = {expr, lb, ub};

                if ~opts.assume_lower_bounds && ~opts.lift_complementarities % Lower bounds on G, H, not already present in MPCC
                    if ~isempty(obj.ind_nonscalar_G)
                        nlp.g.G_lower_bounds = {mpcc.G(obj.ind_nonscalar_G), 0, inf};
                    end
                    if ~isempty(obj.ind_nonscalar_H)
                        nlp.g.H_lower_bounds = {mpcc.H(obj.ind_nonscalar_H), 0, inf};
                    end
                end
                % Get nlpsol plugin
                switch opts.solver
                  case 'ipopt'
                    obj.plugin = nosnoc.solver.plugins.Ipopt();
                  case 'snopt'
                    obj.plugin = nosnoc.solver.plugins.Snopt();
                  case 'worhp'
                    obj.plugin = nosnoc.solver.plugins.Worhp();
                  case 'uno'
                    obj.plugin = nosnoc.solver.plugins.Uno();
                end

                % Construct solver
                obj.plugin.construct_solver(nlp, opts);
                obj.nlp = nlp;
            end
        end

        function met = complementarity_tol_met(obj, stats)
            last_stats = stats.solver_stats(end);
            met = 0;
            if abs(stats.complementarity_stats(end)) < 10 * obj.opts.comp_tol
                met = 1;
            end
        end

        function violation = compute_constraint_violation(obj, w, g)
            nlp = obj.nlp;
            % We get the maximum violation of g: (max(max(lbg - g(w,p), 0), max(g(w,p) - ubg, 0))),
            % and w: (max(max(lbw - w, 0), max(w - ubw, 0)))
            % bounds from the last solve of the relaxed nlp.
            violation = max([nlp.g.mpcc_g().violation; nlp.w.mpcc_w().violation]);
        end

        function out = cat(dim,varargin)
            error('Concatenation not supported')
        end

        function varargout = size(obj,varargin)
            varargout = 1;
        end
        
        function ind = end(obj,k,n)
            ind = 1;
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
            % TODO(@anton) Also perhaps initial complementarity multipliers but need to think about
            %              how exactly to implement that.
            parse(p, index_op(1).Indices{:});

            % data and local variables
            nlp = obj.nlp;
            mpcc = obj.mpcc;
            opts = obj.opts;
            plugin = obj.plugin;
            G_fun = Function('G', {mpcc.x, mpcc.p}, {mpcc.G});
            H_fun = Function('H', {mpcc.x, mpcc.p}, {mpcc.H});
            comp_res_fun = Function('comp_res', {mpcc.x, mpcc.p}, {mmax(mpcc.G.*mpcc.H)});
            f_mpcc_fun = Function('f_mpcc', {mpcc.x, mpcc.p}, {mpcc.f});
            % Update nlp data
            if ~isempty(p.Results.x0)
                nlp.w.mpcc_w().init = p.Results.x0;
            end
            if ~isempty(p.Results.lbx)
                nlp.w.mpcc_w().lb = p.Results.lbx;
            end
            if ~isempty(p.Results.ubx)
                nlp.w.mpcc_w().ub = p.Results.ubx;
            end
            if ~isempty(p.Results.lbg)
                nlp.g.mpcc_g().lb = p.Results.lbg;
            end
            if ~isempty(p.Results.ubg)
                nlp.g.mpcc_g().ub = p.Results.ubg;
            end
            if ~isempty(p.Results.p)
                nlp.p.mpcc_p().val = p.Results.p;
            end
            if ~isempty(p.Results.lam_g0)
                nlp.g.mpcc_g().init_mult = p.Results.lam_g0;
            end
            if ~isempty(p.Results.lam_x0)
                nlp.w.mpcc_w().init_mult = p.Results.lam_x0;
            end

            if ~opts.assume_lower_bounds % Lower bounds on G, H, not already present in MPCC
                lb = nlp.w.mpcc_w().lb;
                lb(obj.ind_map_G) = 0;
                lb(obj.ind_map_H) = 0;
                nlp.w.mpcc_w().lb = lb;
            end
 

            % Initial conditions
            sigma_k = opts.sigma_0;

            % Initialize Stats struct
            stats = struct();
            stats.cpu_time = [];
            stats.wall_time = [];
            stats.cpu_time_total = 0;
            stats.wall_time_total = 0;
            stats.sigma_k = sigma_k;
            stats.homotopy_iterations = [];
            stats.solver_stats = [];
            % TODO(@anton) calculate objective and complementarity residual
            stats.objective = [];
            stats.complementarity_stats = [full(comp_res_fun(nlp.w.mpcc_w().init, nlp.p.mpcc_p().val))];

            % Initialize Results struct
            mpcc_results = struct;
            mpcc_results.W = nlp.w.mpcc_w().init;
            mpcc_results.nlp_results = [];

            % homotopy loop
            complementarity_iter = 1;
            ii = 0;
            last_iter_failed = 0;
            timeout = 0;

            if opts.print_level >= 3
                plugin.print_nlp_iter_header();
            end

            while (abs(complementarity_iter) > opts.comp_tol || last_iter_failed) &&...
                    ii < opts.N_homotopy &&...
                    (sigma_k > opts.sigma_N || ii == 0) &&...
                    ~timeout
                % homotopy parameter update
                last_iter_failed = 0;
                if ii == 0
                    sigma_k = opts.sigma_0;
                else
                    if isequal(opts.homotopy_update_rule,'linear')
                        sigma_k = opts.homotopy_update_slope*sigma_k;
                    elseif isequal(opts.homotopy_update_rule,'superlinear')
                        sigma_k = max(opts.sigma_N,min(opts.homotopy_update_slope*sigma_k,sigma_k^opts.homotopy_update_exponent));
                    else
                        % TODO(@anton) make this not necessary
                        error('For the homotopy_update_rule please select ''linear'' or ''superlinear''.')
                    end
                end
                stats.sigma_k = [stats.sigma_k, sigma_k];
                nlp.p.sigma_p().val = sigma_k;

                if ii ~= 0
                    % TODO(@anton) Lets push on casadi devs to allow for changing of options after construction
                    %              (or maybe do it ourselves) and then remove this hack
                    if obj.opts.timeout_cpu
                        plugin.construct_solver(nlp, opts, opts.timeout_cpu - stats.cpu_time_total);
                    elseif obj.opts.timeout_wall
                        plugin.construct_solver(nlp, opts, opts.timeout_wall - stats.wall_time_total);
                    end % HACK ENDS HERE
                end
                start = cputime;
                tic;
                [nlpsol_stats, nlpsol_results] = nlp.solve();
                cpu_time_iter = cputime - start;
                wall_time_iter = toc;
                nlp.w.init = nlp.w.res;
                
                solver_stats = plugin.cleanup_solver_stats(nlpsol_stats);
                
                stats.solver_stats = [stats.solver_stats, solver_stats];

                mpcc_results.nlp_results = [mpcc_results.nlp_results, nlpsol_results];

                last_iter_failed = plugin.check_iteration_failed(stats);
                timeout = plugin.check_timeout(stats);

                if last_iter_failed
                    obj.print_infeasibility();
                end

                % update timing stats
                stats.cpu_time = [stats.cpu_time,cpu_time_iter];
                stats.cpu_time_total = stats.cpu_time_total + cpu_time_iter;
                stats.wall_time = [stats.wall_time, wall_time_iter];
                stats.wall_time_total = stats.wall_time_total + wall_time_iter;

                if opts.timeout_cpu && (stats.cpu_time_total > opts.timeout_cpu)
                    timeout = 1;
                    last_iter_failed = 1;
                end
                if opts.timeout_wall && (stats.wall_time_total > opts.timeout_wall)
                    timeout = 1;
                    last_iter_failed = 1;
                end
                % update results output.
                mpcc_results.W = [mpcc_results.W,nlp.w.mpcc_w().res]; % all homotopy iterations

                % update complementarity and objective stats
                complementarity_iter = full(comp_res_fun(nlp.w.mpcc_w().res, nlp.p.mpcc_p().val));
                stats.complementarity_stats = [stats.complementarity_stats;complementarity_iter];
                objective = full(f_mpcc_fun(nlp.w.mpcc_w().res, nlp.p.mpcc_p().val));
                stats.objective = [stats.objective, objective];

                % update counter
                ii = ii+1;
                % Verbose
                if opts.print_level >= 3
                    plugin.print_nlp_iter_info(stats)
                end
            end
            
            mpcc_results.f = full(f_mpcc_fun(nlp.w.mpcc_w().res, nlp.p.mpcc_p().val));
            mpcc_results.x = nlp.w.mpcc_w().res;
            mpcc_results.p = nlp.p.mpcc_p().val;
            mpcc_results.lam_x = nlp.w.mpcc_w().mult;
            mpcc_results.g = nlp.g.mpcc_g().eval;
            mpcc_results.lam_g = nlp.g.mpcc_g().mult;
            % TODO(@anton) it may not always be possible to calculate lam_G and lam_H
            %              figure out when this is possible.
            mpcc_results.G = full(G_fun(mpcc_results.x, mpcc_results.p));
            % mpcc_results.lam_G = ;
            mpcc_results.H = full(H_fun(mpcc_results.x, mpcc_results.p));
            % mpcc_results.lam_H = ;
            
            stats.converged = obj.complementarity_tol_met(stats) && ~last_iter_failed && ~timeout;
            stats.success = stats.converged;
            stats.constraint_violation = obj.compute_constraint_violation(mpcc_results.x, mpcc_results.g);

            varargout{1} = mpcc_results;
            obj.stats = stats;
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
