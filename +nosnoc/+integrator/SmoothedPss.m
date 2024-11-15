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

        t_grid
        t_grid_full
        x_res
        x_res_full
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
                fun_opts.allow_free = true;
                f_x_alpha = Function('f_x_alpha', {obj.dcs.alpha}, {obj.dcs.f_x}, fun_opts);
                f_x_smoothed = f_x_alpha(0.5*(1+tanh([obj.model.c{:}]/sigma)));

                rhs_fun = Function('rhs_fun', {model.x, model.u, sigma}, {f_x_smoothed});
                obj.ode_func = @(t, x , u ,sigma) full(rhs_fun(x, u, sigma));
              otherwise
                error("nosnoc: Incorrect model type for smoothed Pss integrator.")
            end
        end
        
        function [t_grid,x_res,t_grid_full,x_res_full] = simulate(obj, plugin, extra_args)
            arguments
                obj nosnoc.integrator.SmoothedPss
                plugin nosnoc.solver.MpccMethod = nosnoc.solver.MpccMethod.SCHOLTES_INEQ
                extra_args.u = []
                extra_args.x0 = [];
            end
            opts = obj.opts;
            solver_opts = obj.solver_opts;
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
            obj.x_res = x0;
            obj.x_res_full = x0;
            obj.t_grid = 0;
            obj.t_grid_full = 0;
            obj.set_x0(x0);
            t_current = 0;

            for ii=1:opts.N_sim
                if ~isempty(extra_args.u) % TODO(@anton) maybe pass as arg to solve
                    u_i = extra_args.u(:,ii);
                else
                    u_i = [];
                end
                switch(obj.solver_opts.matlab_ode_solver)
                  case 'ode45'
                    [t_sim, x_sim] = ode45(@(t, x)  obj.ode_func(t,x,u_i,opts.sigma_smoothing), [t_current, t_current+opts.T], obj.x_curr, solver_opts.matlab_ode_opts);
                  case 'ode23'
                    [t_sim, x_sim] = ode23(@(t, x)  obj.ode_func(t,x,u_i,opts.sigma_smoothing), [t_current, t_current+opts.T], obj.x_curr, solver_opts.matlab_ode_opts);
                  case 'ode113'
                    [t_sim, x_sim] = ode113(@(t, x)  obj.ode_func(t,x,u_i,opts.sigma_smoothing), [t_current, t_current+opts.T], obj.x_curr, solver_opts.matlab_ode_opts);
                  case 'ode78'
                    [t_sim, x_sim] = ode78(@(t, x)  obj.ode_func(t,x,u_i,opts.sigma_smoothing), [t_current, t_current+opts.T], obj.x_curr, solver_opts.matlab_ode_opts);
                  case 'ode89'
                    [t_sim, x_sim] = ode89(@(t, x)  obj.ode_func(t,x,u_i,opts.sigma_smoothing), [t_current, t_current+opts.T], obj.x_curr, solver_opts.matlab_ode_opts);
                  case 'ode15s'
                    [t_sim, x_sim] = ode15s(@(t, x)  obj.ode_func(t,x,u_i,opts.sigma_smoothing), [t_current, t_current+opts.T], obj.x_curr, solver_opts.matlab_ode_opts);
                  case 'ode23s'
                    [t_sim, x_sim] = ode23s(@(t, x)  obj.ode_func(t,x,u_i,opts.sigma_smoothing), [t_current, t_current+opts.T], obj.x_curr, solver_opts.matlab_ode_opts);         
                  case 'ode23t'
                    [t_sim, x_sim] = ode23t(@(t, x)  obj.ode_func(t,x,u_i,opts.sigma_smoothing), [t_current, t_current+opts.T], obj.x_curr, solver_opts.matlab_ode_opts);
                  case 'ode23tb'
                    [t_sim, x_sim] = ode23tb(@(t, x)  obj.ode_func(t,x,u_i,opts.sigma_smoothing), [t_current, t_current+opts.T], obj.x_curr, solver_opts.matlab_ode_opts);
                  case 'ode15i'
                    [t_sim, x_sim] = ode15i(@(t, x)  obj.ode_func(t,x,u_i,opts.sigma_smoothing), [t_current, t_current+opts.T], obj.x_curr, solver_opts.matlab_ode_opts);
                end
                t_current = t_sim(end);
                obj.t_grid = [obj.t_grid, t_sim(end)];
                obj.t_grid_full = [obj.t_grid_full t_sim(2:end)'];
                obj.x_res_full = [obj.x_res_full, x_sim(2:end,:)'];
                obj.x_res = [obj.x_res, x_sim(end,:)'];
                obj.set_x0(x_sim(end,:)');
            end
            x_res = obj.x_res;
            x_res_full = obj.x_res_full;
            t_grid = obj.t_grid;
            t_grid_full = obj.t_grid_full;
        end

        function ret = get(obj, field)
            if ~strcmp(field, 'x')
                error(['nosnoc:' char(field) ' is not a valid field for this integrator.']);
            end
            ret = obj.x_res;
        end

        function ret = get_full(obj, field)
            if ~strcmp(field, 'x')
                error(['nosnoc:' char(field) ' is not a valid field for this integrator.']);
            end
            ret = obj.x_res_full;
        end

        function t_grid = get_time_grid(obj)
            t_grid = obj.t_grid;
        end

        function t_grid_full = get_time_grid_full(obj)
            t_grid_full = obj.t_grid_full;
        end

        function set_x0(obj, x0)
            obj.x_curr = x0;
        end
        
        function clear_history(obj)
            obj.x_res = []; % Clear simulation data.
            obj.x_res_full = []; % Clear simulation data.
            obj.t_grid = []; % Clear simulation data.
            obj.t_grid_full = []; % Clear simulation data.
        end
    end
end
