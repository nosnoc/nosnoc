% This file is part of nosnoc.
classdef Solver < handle & matlab.mixin.indexing.RedefinesParen
    properties
        mpec % Either a struct with the (possibly optional) fields (f, p, w, g, G, H) or a subclass of vdx.problems.Mpcc.
        opts % Options object.
        stats % Struct with stats of the last solve.
        dims % Struct containig some important dimensions
        qpec_problem % Either a struct with the (possibly optional) fields (f, p, w, g, G, H) or a subclass of vdx.problems.Mpcc for storing the parametric QPEC data
        % qpec_solver_options % nosnoc.reg_homotopy.Options or mpecopt.Options for QPEC solver
        qpec_functions % all qpec functions (residuals, first and second order derivatives) (TODO: should live in qpec class later)
        qpec_solver  % nosnoc.reg_homotopy.Solver or nonosc mpccsol is used to solve QPECs
        qpec_initialization % struct storing all qp data. (TODO: should live in qpec class later)
        solver_initialization % initalization of the mpec;
    end

    properties (Access=private)
        iters % solving all iterations of the sqpe method w, lam_g, lam_w, f, inf_primal, inf_dual, d_norm, alpha;
        filter % a struct storing all filter information
        % QPEC functions here??
        %     ind_scalar_G
        %     ind_scalar_H
        %     ind_nonscalar_G
        %     ind_nonscalar_H
        %     ind_map_G % Map from lifted indices to original indices in G.
        %     ind_map_H % Map from lifted indices to original indices in H.
        %
        %     G_fun % Function (nlp.w, nlp.p)|-> mpcc.G
        %     H_fun % Function (nlp.w, nlp.p)|-> mpcc.G
        %     comp_res_fun % Function (nlp.w, nlp.p)|-> mmax(mpcc.G,mpcc.H)
        %     f_mpcc_fun % Function (nlp.w, nlp.p)|-> mpcc.f
        %     w_mpcc_fun % Function (nlp.w)|-> mpcc.w
        %     g_mpcc_fun % Function (nlp.w, nlp.p)|-> mpcc.g
        %
        %     ind_map_g % Map for general constraints containing the corresponding indices in the original MPCC passed in and indices in the NLP
        %     ind_map_w % Map for primal variables containing the corresponding indices in the original MPCC passed in and indices in the NLP
        %     ind_map_p % Map for parameters containing the corresponding indices in the original MPCC passed in and indices in the NLP
    end


    methods (Access=public)

        function obj = Solver(mpec, opts)
            import casadi.*
            t_prepare_qpec = tic;
            obj.opts = opts;
            if isa(mpec, 'vdx.problems.Mpcc')
                obj.mpec = copy(mpec); % copy of mpec
            else
                obj.mpec = mpec;
            end
            % obj.create_mpec_functions();
            % obj.solver_initialization = struct();

            obj.create_qpec_functions();
            obj.create_qpec_solver();
            cpu_time_prepare_qpec = toc(t_prepare_qpec);
            obj.stats.cpu_time_prepare_qpec = cpu_time_prepare_qpec;

            % QPEC specific


        end

        function [solution,stats] = solve(obj, solver_initialization)
            arguments
                obj
                solver_initialization = []
            end
            % Get data
            import casadi.*
            obj.process_solver_initialization(solver_initialization);
            solver_initialization = obj.solver_initialization;
            mpec = obj.mpec;
            opts = obj.opts;
            dims = obj.dims;
            stats = obj.stats;
            obj.iters = struct(); % empty struct for iteration data.

            % main algoritmic loop
            obj.initialize_qpec();
            k = 0; % iteration counter;
            w_k = solver_initialization.x0;
            p0 = solver_initialization.p0;
            %TODO: compute least squares duals
            % TODO@Armin: allow guess for dual variables in solver initalization or do least squares computation
            lam_g_k = zeros(obj.dims.n_g,1);
            lam_w_k = zeros(obj.dims.n_w,1);
            f_k = full(obj.qpec_functions.f_fun(w_k,p0));
            g_k = full(obj.qpec_functions.g_fun(w_k,p0));
            inf_pr_g = max(max(0, max(solver_initialization.lbg - g_k, g_k - solver_initialization.ubg)));
            inf_pr_w = max(max(0, max(solver_initialization.lbx- w_k, w_k - solver_initialization.ubx)));
            inf_pr_k = max(inf_pr_g,inf_pr_w);
            inf_du_k = max(full(obj.qpec_functions.nabla_L_fun(w_k,lam_g_k,lam_w_k,p0)));
            inf_comp_k = max(abs(full(obj.qpec_functions.G_fun(w_k,p0))).*abs(full(obj.qpec_functions.H_fun(w_k,p0))));


            % prealocate dynamic arrays (cells)
            obj.iters.f = {};
            obj.iters.w = {};
            obj.iters.d = {};
            obj.iters.lam_g = {};
            obj.iters.lam_w = {};
            obj.iters.d_norm = {}; % step size norm
            obj.iters.inf_pr = {}; % primal infeasiblity (constarint violation)
            obj.iters.inf_du = {}; % dual infeasiblity (nabla_w L(w,lam_g,lam_w) value)
            obj.iters.inf_comp = {}; % comp infeasiblity
            obj.iters.alpha = {}; % step sizes
            obj.iters.cpu_time_qpec = {}; % step sizes

            % initalize
            obj.iters.k = 0;
            obj.iters.f{end+1} = f_k;
            obj.iters.w{end+1} = w_k;
            obj.iters.d{end+1} = w_k*0;
            obj.iters.lam_g{end+1} = lam_g_k;
            obj.iters.lam_w{end+1} = lam_w_k;
            obj.iters.d_norm{end+1} = 0; % step size norm
            obj.iters.inf_pr{end+1} = inf_pr_k; % primal infeasiblity (constarint violation)
            obj.iters.inf_du{end+1} = inf_du_k; % dual infeasiblity (nabla_w L(w,lam_g,lam_w) value)
            obj.iters.inf_comp{end+1} = inf_comp_k; % dual infeasiblity (nabla_w L(w,lam_g,lam_w) value)
            obj.iters.alpha{end+1} = 1; % % first initla value default 1, not actually used, just a place holder
            obj.iters.cpu_time_qpec{end+1} = 0; % % first initla value default 0, not actually used, just a place holder

            stats.stopping_criterion_fullfiled = false;
            % --------------------------- major SQP itrations -------------------------------------------
            obj.print_iter_header()
            obj.print_iter_stats()
            while (k <= opts.max_iter) && ~stats.stopping_criterion_fullfiled
                % evaluate all functions/create new QPEC
                obj.update_qpec_data();
                % solve QPEC
                t_solve_qpec = tic;
                results_qpec = obj.qpec_solver(obj.qpec_initialization);
                cpu_time_solve_qpec = toc(t_solve_qpec);
                % TODO: exit if QP infeasible - or go to fall back strategy/restoration phase
                d_k = full(results_qpec.x); 
                obj.iters.d{end+1} = d_k; % store currently computed full step
                % if sucess:
                    alpha_k = 1;
                    % do line serach
                    % [alpha_k] = filter_line_serach();
                    % [d_k] = compute_soc();
                    % SOC?
                % else
                    % terminate (no resotration phase yet)
                % end
                
                obj.iters.alpha{end+1} = alpha_k;
                % update iterate
                w_k = w_k + alpha_k*d_k;
                lam_g_k = full(results_qpec.lam_g);
                lam_x_k = full(results_qpec.lam_x);

                % compute all stats:
                f_k = full(obj.qpec_functions.f_fun(w_k,p0));
                g_k = full(obj.qpec_functions.g_fun(w_k,p0));
                inf_pr_g = max(max(0, max(solver_initialization.lbg - g_k, g_k - solver_initialization.ubg)));
                inf_pr_w = max(max(0, max(solver_initialization.lbx- w_k, w_k - solver_initialization.ubx)));
                inf_pr_k = max(inf_pr_g,inf_pr_w);
                inf_du_k = max(full(obj.qpec_functions.nabla_L_fun(w_k,lam_g_k,lam_w_k,p0)));
                inf_comp_k = max(abs(full(obj.qpec_functions.G_fun(w_k,p0))).*abs(full(obj.qpec_functions.H_fun(w_k,p0))));
                d_norm_k = norm(d_k);

                k = k+1;
                % store all stats
                obj.iters.f{end+1} = f_k;
                obj.iters.w{end+1} = w_k;
                obj.iters.lam_g{end+1} = lam_g_k;
                obj.iters.lam_w{end+1} = lam_w_k;
                obj.iters.d_norm{end+1} = d_norm_k; % step size norm
                obj.iters.inf_pr{end+1} = inf_pr_k; % primal infeasiblity (constarint violation)
                obj.iters.inf_du{end+1} = inf_du_k; % dual infeasiblity (nabla_w L(w,lam_g,lam_w) value)
                obj.iters.inf_comp{end+1} = inf_comp_k; % dual infeasiblity (nabla_w L(w,lam_g,lam_w) value)
                obj.iters.alpha{end+1} = alpha_k; % % first initla value default 1, not actually used, just a place holder
                obj.iters.cpu_time_qpec{end+1} = cpu_time_solve_qpec;
                obj.iters.k = k;
                % print details
                obj.print_iter_stats()

                % check convergence criteria


                if d_norm_k <= opts.tol
                    stats.stopping_criterion_fullfiled = true;
                end
            end
            obj.print_iter_line()

            % convert cells to arrays (cells were used because of easier dyn memory allocation)
            obj.iters.f = cell2mat(obj.iters.f);
            obj.iters.w = cell2mat(obj.iters.w);
            obj.iters.lam_w = cell2mat(obj.iters.lam_w);
            obj.iters.lam_g = cell2mat(obj.iters.lam_g);
            obj.iters.d_norm = cell2mat(obj.iters.d_norm);

            solution.x = obj.iters.w(:,end);
            solution.lam_g = obj.iters.lam_g(:,end);
            solution.lam_x = obj.iters.lam_w(:,end);
            solution.f = obj.iters.f(end);
            obj.stats.iters = obj.iters;
            obj.stats.success = true;

        end
        % some overloading
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
    methods(Access=private)

        %  ----------------- algoritmic methods -----------------
        function stopping_criterion_fullfiled = stopping_criterion(obj)
            % TODO
        end

        function process_solver_initialization(obj, solver_initialization)
            import casadi.*
            opts = obj.opts;
            dims = obj.dims;
            qpec_functions = obj.qpec_functions;

            % check does a parameter exist:
            if ~isfield(solver_initialization,"p0")
                solver_initialization.p0 = [];
            end

            G_fun = qpec_functions.G_fun;
            H_fun = qpec_functions.H_fun;
            if opts.initial_comp_all_zero
                G_eval = zeros(dims.n_comp,1);
                H_eval = zeros(dims.n_comp,1);
            else
                G_eval = full(qpec_functions.G_fun(solver_initialization.x0, solver_initialization.p0));
                H_eval = full(qpec_functions.H_fun(solver_initialization.x0, solver_initialization.p0));
            end
            % map_w = dims.map_w;
            % map_g = dims.map_g;
            lbx = solver_initialization.lbx;
            ubx = solver_initialization.ubx;
            lbg = solver_initialization.lbg;
            ubg = solver_initialization.ubg;
            x0 = solver_initialization.x0;
            % TODO@ Armin  CREATE all this residual and helper functions in a seperate method, needed for stopping criteria
            % solver_initialization.lbx = zeros(dims.n_w,1);
            % solver_initialization.lbx(dims.ind_x) = lbx;
            % solver_initialization.lbx(dims.ind_x1) = 0;
            % solver_initialization.lbx(dims.ind_x2) = 0;
            % solver_initialization.ubx = inf(dims.n_w,1);
            % solver_initialization.ubx(dims.ind_x) = ubx;
            % solver_initialization.lbg = zeros(dims.n_g,1);
            % solver_initialization.lbg(dims.ind_g) = lbg;
            % solver_initialization.ubg = zeros(dims.n_g,1);
            % solver_initialization.ubg(dims.ind_g) = ubg;
            % solver_initialization.x0 = zeros(dims.n_w,1);
            % solver_initialization.x0(dims.ind_x) = x0;
            % if opts.lift_complementarities_full
            %     solver_initialization.x0(dims.ind_x1) = G_eval;
            %     solver_initialization.x0(dims.ind_x2) = H_eval;
            % else
            %     solver_initialization.x0(dims.ind_nonscalar_x1) = G_eval(dims.ind_nonscalar_x1);
            %     solver_initialization.x0(dims.ind_nonscalar_x2) = H_eval(dims.ind_nonscalar_x2);
            % end

            % Generate casadi functions for objective and constraint function evaluations
            % Zero order
            % g_eq = g(ind_g_eq)-solver_initialization.lbg(ind_g_eq);                    % g_eq = g - g_lb = 0
            % g_ineq_ub = solver_initialization.ubg(ind_g_ineq_ub)-g(ind_g_ineq_ub);                           % g_ineq_ub = g_ub - g >= 0
            % g_ineq_lb = g(ind_g_ineq_lb)-solver_initialization.lbg(ind_g_ineq_lb);                           % g_ineq_lb = g - g_lb >= 0
            % g_ineq = [solver_initialization.ubg(ind_g_ineq_ub)-g(ind_g_ineq_ub);...
            %     g(ind_g_ineq_lb)-solver_initialization.lbg(ind_g_ineq_lb)];                            % g_ineq = [g_ub; g_lb]
            % n_ineq = size(g_ineq,1);
            % % first-order constraint Jacobians
            % if n_eq > 0
            %     nabla_g_eq = g_eq.jacobian(x);
            % else
            %     nabla_g = [];
            %     nabla_g_eq = [];
            % end
            % if n_ineq > 0
            %     nabla_g_ineq_ub = g_ineq_ub.jacobian(x);
            %     nabla_g_ineq_lb = g_ineq_lb.jacobian(x);
            % else
            %     nabla_g_ineq_ub = [];
            %     nabla_g_ineq_lb = [];
            % end
            % nabla_g_ineq = [nabla_g_ineq_ub; nabla_g_ineq_lb];
            %
            % % Infeasiblity meausre in inf norm
            % % all standard constraints
            % h_eq = max(abs(g_eq));
            % h_ineq_ub = max(min(g_ineq_ub,0));
            % h_ineq_lb = max(min(g_ineq_lb,0));
            % h_ubx = max(min(solver_initialization.ubx-x,0));
            % h_lbx = max(min(x-solver_initialization.lbx,0));
            % % Summary
            % h_std = max([h_eq;h_ineq_ub;h_ineq_lb;h_ubx;h_lbx]);
            % if dims.n_comp > 0
            %     if opts.comp_res_bilinear
            %         h_comp= max(abs(x1).*abs(x2)); % here kappa =0.1 to have reduce of the value by factor of 10
            %                                        % h_comp= sqrt(max(abs(min((x1),(x2))))); % here kappa =0.01 to reduce the value above by factor of 10
            %     else
            %         h_comp= max(min(abs(x1),abs(x2))); %
            %         % h_comp = max(abs(min(x1,x2)));
            %     end
            % else
            %     h_comp  = 0;
            % end
        end

        % problem set up methods
        function create_qpec_functions(obj)
            import casadi.*
            if isa(obj.mpec, 'vdx.problems.Mpcc')
                % get all symbolic expressions
                w = obj.mpec.w.sym;
                p = obj.mpec.p.sym;
                f = obj.mpec.f;
                g = obj.mpec.g.sym;
                G = obj.mpec.G.sym;
                H = obj.mpec.H.sym;
            else
                % get all symbolic expressions
                try
                    w = obj.mpec.w;
                catch
                    w = obj.mpec.x;
                end
                p = obj.mpec.p;
                f = obj.mpec.f;
                g = obj.mpec.g;
                G = obj.mpec.G;
                H = obj.mpec.H;
            end

            % Get some dimensions;
            obj.dims.n_g = length(g);
            obj.dims.n_w = length(w);
            obj.dims.n_comp = length(G);
            obj.dims.n_p = length(p);

            % Define functions for derivatives

            lam_g = SX.sym('lam_g', obj.dims.n_g); % Lagrange multiplier for all general constraints
            lam_w = SX.sym('lam_g', obj.dims.n_w); % Lagrange multiplier for all box constraints
            if  obj.dims.n_g > 0
                L = f+lam_g'*g+lam_w'*w;
            else
                L = f+lam_w'*w;
            end

            % objective
            obj.qpec_functions.f_fun = Function('f', {w, p}, {f});
            obj.qpec_functions.nabla_f_fun = Function('nabla_f', {w, p}, {gradient(f, w)});
            obj.qpec_functions.nabla_L_fun = Function('nabla_L', {w, lam_g, lam_w, p}, {jacobian(L, w)});
            obj.qpec_functions.hess_L_fun = Function('hess_L', {w, lam_g, p}, {hessian(L, w)});
            % regular constraints
            obj.qpec_functions.g_fun = Function('g', {w, p}, {g});
            obj.qpec_functions.nabla_g_fun = Function('nabla_g', {w, p}, {jacobian(g, w)});

            % complementarity constraints;
            obj.qpec_functions.G_fun = Function('G', {w, p}, {G});
            obj.qpec_functions.H_fun = Function('H', {w, p}, {H});
            obj.qpec_functions.nabla_G_fun = Function('nabla_G', {w, p}, {jacobian(G, w)});
            obj.qpec_functions.nabla_H_fun = Function('nabla_H', {w, p}, {jacobian(H, w)});
            % obj.qpec_functions.w = w;
        end

        function create_qpec_solver(obj)
            % TODO@Armin: MX or SX better here, especially for matrices? Which type?
            import casadi.*
            % general problem parameters:
            if isa(obj.mpec, 'vdx.problems.Mpcc')
                p = obj.mpec.p.sym;
            else
                p = obj.mpec.p;
            end


            % Parameters for QPEC vectors and matrices
            d = SX.sym('d', obj.dims.n_w, 1); % dof in QPEC
            % Objective data parameters
            nabla_f_p = SX.sym('nabla_f', obj.dims.n_w, 1);
            Q_p = SX.sym('Q_p', obj.qpec_functions.hess_L_fun.sparsity_out(0));
            % Regular constraint data parameters
            g_p = SX.sym('g_p', obj.dims.n_g, 1);
            nabla_g_p = SX.sym('nabla_g_p', obj.qpec_functions.nabla_g_fun.sparsity_out(0));
            % Complementarity constraints data parameters
            G_p = SX.sym('G_p', obj.dims.n_comp, 1);
            nabla_G_p = SX.sym('nabla_G_p', obj.qpec_functions.nabla_G_fun.sparsity_out(0));
            H_p = SX.sym('H_p', obj.dims.n_comp, 1);
            nabla_H_p = SX.sym('nabla_H_p', obj.qpec_functions.nabla_H_fun.sparsity_out(0)); %
            % Algorithmic parameters
            rho_d = SX.sym('rho_d', 1); % penalty parameter on step length d: rho_h\|d\|_2^2
            % all problem data
            p_qpec = vertcat(Q_p(:),nabla_f_p,nabla_g_p(:),nabla_G_p(:),nabla_H_p(:),g_p,G_p,H_p,rho_d);

            % QPEC problem expressions
            f_k = 0.5*d'*Q_p*d+ nabla_f_p'*d+rho_d*d'*d;
            g_k = g_p + nabla_g_p*d;
            G_k = G_p + nabla_G_p*d;
            H_k = H_p + nabla_H_p*d;
            % QPEC structure for mpccsol
            obj.qpec_problem = struct();
            obj.qpec_problem.x = d;
            obj.qpec_problem.f = f_k;
            obj.qpec_problem.g = g_k;
            obj.qpec_problem.G = G_k;
            obj.qpec_problem.H = H_k;
            obj.qpec_problem.p = vertcat(p_qpec,p); % all problem parameters

            if isa(obj.opts.qpec_solver_options,'nosnoc.reg_homotopy.Options')
                obj.qpec_solver = nosnoc.solver.mpccsol('qpec_solver', 'reg_homotopy', obj.qpec_problem, obj.opts.qpec_solver_options);
            else
                obj.qpec_solver = nosnoc.solver.mpccsol('qpec_solver', 'mpecopt', obj.qpec_problem, obj.opts.qpec_solver_options);
            end
        end


        function initialize_qpec(obj)
            % do basic initalization of qpec;
            obj.qpec_initialization = struct();
            % need to read out of solver initalization:
            d0 = zeros(obj.dims.n_w,1); %
            % Bounds on general constriants
            % (no shift, as constraints read as lb <= g_k + nabla_g'*w <= ub)
            obj.qpec_initialization.lbg = obj.solver_initialization.lbg;
            obj.qpec_initialization.ubg = obj.solver_initialization.ubg;
            obj.qpec_initialization.x0 = d0; % initial guess for qpecs solver
            % default bounds on w
            obj.qpec_initialization.lbx = obj.solver_initialization.lbx - obj.solver_initialization.lbx;
            obj.qpec_initialization.ubx = obj.solver_initialization.ubx - obj.solver_initialization.lbx;
        end


        function update_qpec_data(obj)
            % function update_qpec(obj)
            % This method does all funciton evaluations needed for constructing a new QPEC
            % Read most up to date data
            p_val = obj.solver_initialization.p0;
            w_k = obj.iters.w{end}; % @TODO: This should be stored somewhere else, e.g., iters and read the last one

            if obj.opts.discard_constraints_in_hessian
                lam_g_k = zeros(obj.dims.n_g);
            else
                lam_g_k = obj.iters.lam_g{end};
            end

            % Objective
            nabla_f_k = full(obj.qpec_functions.nabla_f_fun(w_k, p_val));
            Q_k = full(obj.qpec_functions.hess_L_fun(w_k, lam_g_k, p_val));
            % Regular constraints
            g_k = full(obj.qpec_functions.g_fun(w_k, p_val));
            nabla_g_k = full(obj.qpec_functions.nabla_g_fun(w_k, p_val));
            % Comp G
            G_k = full(obj.qpec_functions.G_fun(w_k, p_val));
            nabla_G_k = full(obj.qpec_functions.nabla_G_fun(w_k, p_val));
            % Comp H
            H_k = full(obj.qpec_functions.H_fun(w_k, p_val));
            nabla_H_k = full(obj.qpec_functions.nabla_H_fun(w_k, p_val));
            % all problem data
            p_qpec_k = vertcat(Q_k(:),nabla_f_k,nabla_g_k(:),nabla_G_k(:),nabla_H_k(:),g_k,G_k,H_k,obj.opts.rho_d,p_val);
            % problem parameters;
            obj.qpec_initialization.p = p_qpec_k;

            % Shift bounds  as lb <= w_k+d <= ub reads as lb - w_k <= d <= ub - w_k
            obj.qpec_initialization.lbx = obj.solver_initialization.lbx - w_k;
            obj.qpec_initialization.ubx = obj.solver_initialization.ubx - w_k;
        end

        function alpha = filter_line_serach(obj,d_qpec)
            % output accept_trail_step, filter, filter_message] = ...
            % check_is_acceptable_to_filter(x_k, x_k_trail, rho_TR_k_l, filter, p_val, settings)
            % Biegler 2010  - Algorithm 5.4, p. 120-121;
            % TODO described filter.filter_entries

            % Filter check
            p_val = obj.solver_initialization.p0;
            w_k = obj.iters.w{end}; 
            d_k = obj.iters.d{end}; 
            w_k_trail = w_k+alpha_k*d_k;
            % Current iterate - TODO @ Read this from QP DATA as already evaluated
            % f_k = full(filter.f_fun(x_k, p_val));
            
            % h_std_k = full(filter.h_std_fun(x_k, p_val));
            % h_comp_k = full(filter.h_comp_fun(x_k, p_val));
            % h_total_k = max(h_std_k,h_comp_k);
            % Trail point  - TODO: e
            f_k_trail = full(obj.qpec_functions.f_fun(w_k_trail, p_k));
            % h_std_k_trail = full(filter.h_std_fun(x_k_trail, p_val));
            % h_comp_k_trail = full(filter.h_comp_fun(x_k_trail, p_val));
            % h_total_k_trail = max(h_std_k_trail,h_comp_k_trail);
            % f_actual_reduction = f_k - f_k_trail;

            % Directional derivative
            % m_k_l  = nabla_f(x_k)^T d_k;
            % m_k_l = -rho_TR_k_l*filter.f_predicted_reduction;
            m_k_l = nabla_f_k'*d_k;
            % gamma_f = 1e-3;
            % gamma_h = 1e-3;
            % f-swtiching condition constantsl;
            % h_max = 1e4*max(1,h(x_0));
            % h_min = 1e-4*min(1,h(x_0));
            % if rho_TR_k_l too small  - restoration phase,
            % 4 (c) -  Check if acceptabel to filter:
            % acceptable_to_filter = false;
            % for ii = 1:length(filter.filter_entries)
            %     if f_k_trail <= filter.filter_entries{ii}(1) || h_std_k_trail <= filter.filter_entries{ii}(4)|| h_comp_k_trail <= filter.filter_entries{ii}(4)
            %         acceptable_to_filter  = true;
            %     else
            %         acceptable_to_filter  = false;
            %         break;
            %     end
            % end

            % if acceptable_to_filter
            %     % Check if sufficient objective decerase w.r.t to current iterate
            %     feasible_enough = h_total_k < settings.filter_h_min; % Check A 5.8. and Eq. (19) in Waecher ipopt paper;
            %     f_type_switching_codnitions_holds = feasible_enough &&  m_k_l < 0 && (-m_k_l)^settings.filter_s_f*(rho_TR_k_l)^(1-settings.filter_s_f) > settings.filter_delta*(h_total_k)^settings.filter_s_h;
            %     if f_type_switching_codnitions_holds
            %         % Case I (f-type): does switching condition hold: TRUE
            %         % & does armijo
            %         armijo_condition_holds = f_k_trail <= f_k + settings.filter_eta_f*m_k_l;
            %         % armijo_condition_holds =  f_k - f_k_trail >= -settings.filter_eta_f*m_k_l;
            %         if armijo_condition_holds
            %             accept_trail_step = true;
            %             filter_message = 'f';
            %         else
            %             accept_trail_step = false;
            %         end
            %     else
            %         % Case II (h-type);
            %         % IF sufficient decerase w.r.t CURRENT iterate, i.e. slanting filter condition
            %         slanting_filter_condition_holds = (h_total_k_trail <=(1-settings.filter_gamma_h)*h_total_k) || (f_k_trail <=f_k-settings.filter_gamma_f*h_total_k_trail);
            %         if slanting_filter_condition_holds
            %             accept_trail_step = true;
            %             filter_message = 'h';
            %             % augment filter
            %             filter.Filter_entries = [filter.Filter_entries, [f_k_trail,h_std_k_trail,h_comp_k_trail,h_total_k_trail]];
            %         else
            %             accept_trail_step = false;
            %         end
            %     end
            % else
            %     accept_trail_step = false;
            % end
            % % 4(e)
            % % reject step, reduce TR and compute new step (check the filte only if nlp resolved)
            % if ~accept_trail_step
            %     % filter.Filter_entries = [filter.Filter_entries,[f_k_trail,h_std_k_trail,h_comp_k_trail,h_total_k_trail]]; % augment filter     % will add lot of entries to the filter;
            %     filter_message = 'reject';
            % end
        end

        %  ----------------- verbose methods -----------------
        function print_iter_header(obj)
            line_str = repmat('-',1,97);
            fprintf(['\n' line_str])
            fprintf('\n|%-7s|%-15s|%-10s|%-10s|%-10s|%-10s|%-10s|%-7s|%-7s|\n', ...
                'iter', 'objective', 'inf_pr', 'inf_du' , 'inf_comp', '||d||', 'alpha' ,'Time (s)','Ls.');
            fprintf([line_str '\n'])

        end

        function print_iter_line(obj)
            line_str = repmat('-',1,97);
            line_str = [line_str '\n'];
            fprintf(line_str)
        end

        function print_iter_stats(obj)
            fprintf('|%-7d|%-15.2e|%-10.2e|%-10.2e|%-10.2e|%-10.2e|%-10.2e|%-7.2e|\n', ...
                obj.iters.k, obj.iters.f{end}, obj.iters.inf_pr{end}, obj.iters.inf_du{end}, obj.iters.inf_comp{end}, obj.iters.d_norm{end},  obj.iters.alpha{end}, obj.iters.cpu_time_qpec{end});
        end

        function [] = print_iter_summary(f,h_std,h_comp,solver_message,multiplier_based_stationarity,b_stationariry,n_biactive,f_lpec,rho_TR)

            % print_iter_line()
            line_str = repmat('-',1,130);
            line_str = ['\n' line_str '\n' ];
            fprintf(line_str)
            fprintf('MPECopt:\t %s\n',solver_message);
            fprintf('Objective....................:\t %2.6e\n',f);
            fprintf('Std. constraint violation....:\t %2.6e\n',h_std);
            fprintf('Complementarity residual.....:\t %2.6e\n',h_comp);
            % fprintf('Solver message...............:\t %s\n',solver_message)
            fprintf('Mult. based stationarity.....:\t %s\n',multiplier_based_stationarity);
            fprintf('B stationarity...............:\t %s\n',mat2str(b_stationariry));
            fprintf('Biactive constraints.........:\t %d\n',n_biactive);
            fprintf('nabla_f(x)^T d...............:\t %d\n',f_lpec);
            fprintf('Final rho_TR.................:\t %d\n',rho_TR);
            print_iter_line()
            fprintf('\n');
        end

    end



    % functional methods
    methods(Access=protected)
        function varargout = parenReference(obj, index_op)
            import casadi.*;
            p = inputParser;
            addParameter(p, 'x0', []);
            addParameter(p, 'y0', []);
            addParameter(p, 'lbx', []);
            addParameter(p, 'ubx', []);
            addParameter(p, 'lbg', []);
            addParameter(p, 'ubg', []);
            addParameter(p, 'p', []);
            addParameter(p, 'lam_g0', []);
            addParameter(p, 'lam_x0', []);
            parse(p, index_op(1).Indices{:});

            dims = obj.dims;

            if ~isempty(p.Results.x0)
                obj.solver_initialization.x0 = p.Results.x0;
            else
                obj.solver_initialization.x0 = zeros(dims.n_w,1);
            end
            if ~isempty(p.Results.lbx)
                obj.solver_initialization.lbx = p.Results.lbx;
            else
                obj.solver_initialization.lbx = -inf(dims.n_w,1);
            end
            if ~isempty(p.Results.ubx)
                obj.solver_initialization.ubx = p.Results.ubx;
            else
                obj.solver_initialization.ubx = inf(dims.n_w,1);
            end
            if ~isempty(p.Results.lbg)
                obj.solver_initialization.lbg = p.Results.lbg;
            else
                obj.solver_initialization.lbg = zeros(dims.n_g,1);
            end
            if ~isempty(p.Results.ubg)
                obj.solver_initialization.ubg = p.Results.ubg;
            else
                obj.solver_initialization.ubg = zeros(dims.n_g,1);
            end
            if ~isempty(p.Results.p)
                obj.solver_initialization.p0 = p.Results.p;
            else
                try
                    obj.solver_initialization.p0 = zeros(dims.n_p,1);
                catch
                    obj.solver_initialization.p = zeros(dims.n_p,1);

                end
            end
            if ~isempty(p.Results.lam_g0)
                obj.solver_initialization.lam_g0 = p.Results.lam_g0;
            else
                obj.solver_initialization.lam_g0 = zeros(dims.n_g, 1);
            end
            if ~isempty(p.Results.lam_x0)
                obj.solver_initialization.lam_x0 = p.Results.lam_x0;
            else
                obj.solver_initialization.lam_x0 = zeros(dims.n_w, 1);
            end
            if ~isempty(p.Results.y0)
                obj.solver_initialization.y0 = p.Results.y0;
            end

            [solution,stats] = obj.solve(obj.solver_initialization);
            varargout{1} = solution;
            obj.stats = stats;
        end

        function obj = parenAssign(obj,index_op,varargin)
            % nosnoc.error('invalid', 'Invalid operation');
            % TODO: there is no nosnoc in mpecopt - adjust errors messages
            error('mpecopt: Invalid operation.')
        end

        function obj = parenDelete(obj,index_op)
            % nosnoc.error('invalid', 'Invalid operation')
            % TODO: there is no nosnoc in mpecopt - adjust errors messages
            error('mpecopt: Invalid operation.')
        end

        function n = parenListLength(obj,index_op,ctx)
            n = 1;
        end
    end
end