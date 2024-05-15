classdef solver < handle
    properties % TODO separate these by Get/Set access
        model
        opts
        solver_opts

        dcs
        discrete_time_problem
        stats
    end

    methods
        function obj = solver(model, opts, solver_opts)
            obj.model = model;
            obj.opts = opts;
            obj.solver_opts = solver_opts;

            % Always process model and options
            opts.preprocess();
            model.verify_and_backfill(opts);

            % Run pipeline
            switch class(model)
              case "nosnoc.model.pss"
                if opts.dcs_mode == DcsMode.Stewart
                    obj.dcs = nosnoc.dcs.stewart(model);
                    obj.dcs.generate_variables(opts);
                    obj.dcs.generate_equations(opts);
                    obj.discrete_time_problem = nosnoc.discrete_problem.stewart(obj.dcs, opts);
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

        function x = getX(obj)
            opts = obj.opts;
            if opts.right_boundary_point_explicit
                x = obj.discrete_time_problem.w.x(:,:,obj.opts.n_s).res;
            else
                x = [obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).res,...
                    obj.discrete_time_problem.w.x(1:opts.N_stages,1:opts.N_finite_elements(1),obj.opts.n_s+1).res];
            end
        end

        function u = getU(obj)
            u = obj.discrete_time_problem.w.u(:).res;
        end

        function t_grid = getTimeGrid(obj)
            h = obj.discrete_time_problem.w.h(:,:).res;
            t_grid = cumsum([0, h]);
        end

        function t_grid = getControlGrid(obj)
            t_grid = [0];
            for ii=1:obj.opts.N_stages
                h_sum = sum(obj.discrete_time_problem.w.h(ii,:).res);
                t_grid = [t_grid, t_grid(end)+h_sum];
            end
        end

        function results = getResults(obj)
            results = obj.discrete_time_problem.get_results_struct();
        end
    end
end
