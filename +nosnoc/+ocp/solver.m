classdef solver < handle
    properties % TODO separate these by Get/Set access
        model
        opts
        solver_opts

        dcs
        discrete_problem
        stats
    end

    methods
        function obj = solver(model, opts, solver_opts)
            obj.model = model;
            obj.opts = opts;
            obj.solver_opts = solver_opts;

            switch class(model)
              case "nosnoc.model.pss"
                if opts.dcs_mode == DcsMode.Stewart
                    obj.dcs = nosnoc.dcs.stewart(model);
                    obj.discrete_problem = nosnoc.discrete_problem.stewart(obj.dcs, opts);
                    obj.discrete_problem.populate_problem();
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
            obj.discrete_problem.create_solver(obj.solver_options, plugin);

            obj.stats = obj.discrete_problem.solve();
        end

        function x = getX(obj)
            x = obj.discrete_problem.w.u(:,:,obj.opts.n_s).res;
        end

        function u = getU(obj)
            u = obj.discrete_problem.w.u(:).res;
        end
    end
end
