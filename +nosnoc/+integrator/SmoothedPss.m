classdef SmoothedPss < handle
    properties
        model
        opts
        solver_opts

        dcs
        ode_func
    end

    properties(Access=private)
        x_curr

        x_all
    end

    methods
        function obj = SmoothedPss(model, opts, solver_opts)
            import casadi.*
            obj.model = model;
            obj.opts = opts;
            obj.solver_opts = solver_opts;

            opts.preprocess();
            model.verify_and_backfill(opts);

            if class(model) == "nosnoc.model.Cls" && opts.time_freezing
                model = nosnoc.time_freezing.reformulation(model, opts);
                model.verify_and_backfill(opts);
                obj.model = model;
            end

            % Run pipeline
            switch class(model)
              case "nosnoc.model.Pss"
                obj.dcs = nosnoc.dcs.Heaviside(model);
                obj.dcs.generate_variables(opts);
                obj.dcs.generate_equations(opts);

                sigma = SX.sym('sigma');
                f_x_alpha = Function('f_x_alpha', {obj.dcs.alpha}, {obj.dcs.f_x});
                f_x_smoothed = f_x_alpha(tanh(sigma*[obj.model.c{:}]));

                rhs_fun = Function('rhs_fun', {model.x, model.u, sigma}, {f_x_smoothed});
                obj.ode_func = @(t, x , u ,sigma) full(rhs_fun(x, u, sigma));
              otherwise
                error("nosnoc: Incorrect model type for smoothed Pss integrator.")
            end
        end
        
        function [t_grid,x_res] = simulate(obj, plugin, extra_args)
            arguments
                obj nosnoc.integrator.SmoothedPss
                plugin nosnoc.solver.MpccMethod = nosnoc.solver.MpccMethod.SCHOLTES_INEQ
                extra_args.u = []
                extra_args.x0 = [];
            end
            opts = obj.opts;
            % TODO(@anton) validators here.
            if ~isempty(extra_args.u) && any(size(extra_args.u) ~= [obj.model.dims.n_u opts.N_sim])
                error("nosnoc: wrong dimensions passed for controls to integrator.")
            end
            if ~isempty(extra_args.x0)
                if  any(size(extra_args.x0) ~= size(obj.model.x0))
                    error("nosnoc: wrong dimensions passed for controls to integrator.")
                end
                x0 = extra_args.x0;
            else
                x0 = obj.model.x0;
            end
            x_res = x0;
            t_grid = 0;
            obj.set_x0(x0);
            t_current = 0;

            for ii=1:opts.N_sim
                if ~isempty(extra_args.u) % TODO(@anton) maybe pass as arg to solve
                    u_i = extra_args.u(:,ii);
                else
                    u_i = [];
                end
                [t_current, x_sim] = ode23s(@(t, x)  obj.ode_func(t,x,u_i,opts.sigma_smoothing), [t_current, t_current+opts.T], obj.x_curr); % todo add options

                t_grid = [t_grid t_current(2:end)'];
                obj.x_all = [obj.x_all, x_sim'];
                obj.set_x0(x_sim);
            end
            x_res = obj.x_all;
        end

        function set_x0(obj, x0)
            obj.x_curr = x0;
        end
        
        function clear_history(obj)
            obj.x_all = []; % Clear simulation data.
        end
    end
end
