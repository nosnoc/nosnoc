classdef Integrator < handle

    properties
        plugin
        model
        opts
        integrator_opts
    end

    methods
        function obj = Integrator(model, opts, integrator_opts)
            obj.model = model;
            obj.opts = opts;
            obj.integrator_opts = integrator_opts;
            integrator_opts.preprocess();
            opts.preprocess();
            model.verify_and_backfill(opts);

            switch integrator_opts.integrator_plugin
              case IntegratorType.FESD
                obj.plugin = nosnoc.integrator.FESD(model, opts, integrator_opts);
              case IntegratorType.SMOOTHED_PSS
                obj.plugin = nosnoc.integrator.SmoothedPss(model, opts, integrator_opts);
            end
        end

        function [t_grid,x_res,t_grid_full,x_res_full] = simulate(obj, extra_args)
            arguments
                obj nosnoc.Integrator
                extra_args.u = []
                extra_args.x0 = [];
            end
            [t_grid,x_res,t_grid_full,x_res_full] = obj.plugin.simulate(u=extra_args.u, x0=extra_args.x0);
            
        end

        function ret = get(obj, field)
            ret = obj.plugin.get(field);
        end

        function ret = get_full(obj, field)
            ret = obj.plugin.get_full(field);
        end

        function t_grid = get_time_grid(obj)
            t_grid = obj.plugin.get_time_grid();
        end

        function t_grid_full = get_time_grid_full(obj)
            t_grid_full = obj.plugin.get_time_grid_full();
        end

        function set(obj, varname, field, indices, value)
            obj.plugin.set(varname, field, indices, value);
        end

        function set_x0(obj, x0)
            obj.plugin.set_x0(x0);
        end

        function set_param(obj, param, value)
            obj.plugin.set_param(param, value);
        end
    end
end
