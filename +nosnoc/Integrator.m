classdef Integrator < handle
    properties % TODO separate these by Get/Set access
        model
        opts
        solver_opts

        dcs
        discrete_time_problem
        stats
    end

    properties(Access=private)
        w_all
    end

    methods
        function obj = Integrator(model, opts, solver_opts)
            obj.model = model;
            obj.opts = opts;
            obj.solver_opts = solver_opts;

            % Always process model and options
            % for integrator also take the extra step of re-calculating N_stages/N_finite_elements
            opts.preprocess();
            if opts.N_stages > 1
                warning("Integrator created with more than 1 control stage. Converting this to finite elements.")
                N_fe = sum(opts.N_finite_elements)
                opts.N_finite_elements = N_fe;
                opts.N_stages = 1;
            end
            model.verify_and_backfill(opts);

            % Run pipeline
            switch class(model)
              case "nosnoc.model.Pss"
                if opts.dcs_mode == DcsMode.Stewart
                    obj.dcs = nosnoc.dcs.Stewart(model);
                    obj.dcs.generate_variables(opts);
                    obj.dcs.generate_equations(opts);
                    obj.discrete_time_problem = nosnoc.discrete_time_problem.Stewart(obj.dcs, opts);
                    obj.discrete_time_problem.populate_problem();
                elseif opts.dcs_mode == DcsMode.Step % TODO: RENAME
                    error("not implemented")
                else
                    error("PSS models can only be reformulated using the Stewart or Heaviside Step reformulations.")
                end
              case "nosnoc.model.heaviside"
                error("not implemented")
              case "nosnoc.model.cls"
                error("not implemented")
              case "nosnoc.model.pds"
                error("not implemented")
            end
        end

        function solve(obj, plugin)
            if ~exist('plugin', 'var')
                plugin = 'scholtes_ineq';
            end
            obj.discrete_time_problem.create_solver(obj.solver_opts, plugin);

            obj.stats = obj.discrete_time_problem.solve();
        end

        function [t_grid,x_res] = simulate(obj)
            if ~exist('plugin', 'var')
                plugin = 'scholtes_ineq';
            end
            opts = obj.opts;
            x_res = obj.model.x0;
            t_grid = 0;
            obj.set_x0(obj.model.x0);
            obj.w_all = []; % Clear simulation data.

            for ii=1:opts.N_sim
                obj.solve(plugin);
                obj.w_all = [obj.w_all,obj.discrete_time_problem.w.res];
                x_step = obj.get_x();
                x_res = [x_res, x_step(:,2:end)];
                if opts.use_fesd
                    h = obj.discrete_time_problem.w.h(:,:).res;
                else
                    h = ones(1,opts.N_finite_elements) * obj.discrete_time_problem.p.T().val/opts.N_finite_elements;
                end
                t_grid = [t_grid, t_grid(end) + cumsum(h)];
                obj.set_x0(x_step(:,end));
            end
        end

        function ret = get(obj, field)
            if isempty(obj.w_all)
                ret = [];
                return
            end

            opts = obj.opts;
            try
                var = obj.discrete_time_problem.w.(field);
            catch
                error([char(field) ' is not a valid field for this integrator.']);
            end

            obj.discrete_time_problem.w.res = obj.w_all(:,1);
            if var.depth == 3
                if opts.right_boundary_point_explicit
                    ret = var(:,:,obj.opts.n_s).res;
                else
                    ret = [var(0,0,obj.opts.n_s).res,...
                        var(1:opts.N_stages,1:opts.N_finite_elements(1),obj.opts.n_s+1).res];
                end
            else
                indexing(1:var.depth) = {':'};
                ret = var(indexing{:}).res;
            end

            for ii=2:size(obj.w_all, 2)
                obj.discrete_time_problem.w.res = obj.w_all(:,ii);
                if var.depth == 3
                    if opts.right_boundary_point_explicit
                        ret = [ret, var(1,:,obj.opts.n_s).res];
                    else
                        ret = [ret, var(1,1:opts.N_finite_elements(1),obj.opts.n_s+1).res];
                    end
                else
                    indexing(2:var.depth) = {':'};
                    if var.depth
                        indexing(1) = {1};
                    end
                    ret = [ret,var(indexing{:}).res];
                end 
            end
        end

        function ret = get_full(obj, field)
            if isempty(obj.w_all)
                ret = [];
                return
            end

            opts = obj.opts;
            try
                var = obj.discrete_time_problem.w.(field);
            catch
                error([char(field) ' is not a valid field for this integrator.']);
            end

            obj.discrete_time_problem.w.res = obj.w_all(:,1);
            indexing(2:var.depth) = {':'};
            if var.depth
                indexing(1) = {1};
            end
            ret = var(indexing{:}).res;

            for ii=2:size(obj.w_all, 2)
                obj.discrete_time_problem.w.res = obj.w_all(:,ii);
                indexing(2:var.depth) = {':'};
                if var.depth
                    indexing(1) = {1};
                end
                ret = [ret,var(indexing{:}).res];
            end
        end

        function t_grid = get_time_grid(obj)
            opts = obj.opts;
            if opts.use_fesd
                h = obj.get('h');
            else
                h = [];
                for ii=1:size(obj.w_all,2)
                    h = [h,ones(1,opts.N_finite_elements) * obj.discrete_time_problem.p.T().val/opts.N_finite_elements];
                end
            end
            t_grid = cumsum([0, h]);
        end

        function x = get_x(obj)
            opts = obj.opts;
            if opts.right_boundary_point_explicit
                x = obj.discrete_time_problem.w.x(:,:,obj.opts.n_s).res;
            else
                x = [obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).res,...
                    obj.discrete_time_problem.w.x(1:opts.N_stages,1:opts.N_finite_elements(1),obj.opts.n_s+1).res];
            end
        end

        function x = get_xend(obj)
            if opts.right_boundary_point_explicit
                x = obj.discrete_time_problem.w.x(1,opts.N_finite_elements,obj.opts.n_s).res;
            else
                x = obj.discrete_time_problem.w.x(1,opts.N_finite_elements,obj.opts.n_s+1).res;
            end
        end

        function set_x0(obj, x0)
            obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).init = x0;
            obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).lb = x0;
            obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).ub = x0;
        end
    end
end
