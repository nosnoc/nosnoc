classdef Solver < handle
    properties % TODO separate these by Get/Set access.
        model % Model for the current OCP, after any pre-processing.
        opts % Nosnoc/FESD discretization options.
        solver_opts % Options passed to MPCC solver.

        dcs % Dynamic Complementarity syste into which the model is processed.
        discrete_time_problem % Discretized optimal control problem.
        stats % Stats, containing both solver and reformulation information.

        active_set % Current active set. TODO(@anton) do we need to store this? Should it be updated after solve? etc.
                   % TODO(@anton) do we want to default this to x0 initialization?
    end

    methods
        function obj = Solver(model, opts, solver_opts)
            obj.model = model;
            obj.opts = opts;
            obj.solver_opts = solver_opts;

            % Always process model and options
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
                if opts.dcs_mode == DcsMode.Stewart
                    obj.dcs = nosnoc.dcs.Stewart(model);
                    obj.dcs.generate_variables(opts);
                    obj.dcs.generate_equations(opts);
                    obj.discrete_time_problem = nosnoc.discrete_time_problem.Stewart(obj.dcs, opts);
                    obj.discrete_time_problem.populate_problem();
                elseif opts.dcs_mode == DcsMode.Heaviside
                    obj.dcs = nosnoc.dcs.Heaviside(model);
                    obj.dcs.generate_variables(opts);
                    obj.dcs.generate_equations(opts);
                    obj.discrete_time_problem = nosnoc.discrete_time_problem.Heaviside(obj.dcs, opts);
                    obj.discrete_time_problem.populate_problem();
                else
                    nosnoc.error('wrong_dcs_mode', "PSS models can only be reformulated using the Stewart or Heaviside Step reformulations.")
                end
              case "nosnoc.model.Heaviside"
                obj.dcs = nosnoc.dcs.Heaviside(model);
                obj.dcs.generate_variables(opts);
                obj.dcs.generate_equations(opts);
                obj.discrete_time_problem = nosnoc.discrete_time_problem.Heaviside(obj.dcs, opts);
                obj.discrete_time_problem.populate_problem();
              case "nosnoc.model.Cls"
                if ~opts.use_fesd
                    nosnoc.error('fesd_j_without_fesd',"The FESD-J reformulation only makes sense with use_fesd=true.")
                end
                obj.dcs = nosnoc.dcs.Cls(model);
                obj.dcs.generate_variables(opts);
                obj.dcs.generate_equations(opts);
                obj.discrete_time_problem = nosnoc.discrete_time_problem.Cls(obj.dcs, opts);
                obj.discrete_time_problem.populate_problem();
              case "nosnoc.model.Pds"
                if ~opts.right_boundary_point_explicit
                    nosnoc.error('pds_rbp_not_one', "You are using an rk scheme with its right boundary point (c_n) not equal to one. Please choose another scheme e.g. RADAU_IIA.")
                end
                if opts.rk_representation == RKRepresentation.differential
                    nosnoc.error('pds_differential',"Differential representation without lifting is unsupported for gradient complementarity systems. Use integral or lifted differential representation.")
                end
                obj.dcs = nosnoc.dcs.Gcs(model);
                obj.dcs.generate_variables(opts);
                obj.dcs.generate_equations(opts);
                obj.discrete_time_problem = nosnoc.discrete_time_problem.Gcs(obj.dcs, opts);
                obj.discrete_time_problem.populate_problem();
              case "nosnoc.model.PDSObjects"
                if ~opts.right_boundary_point_explicit
                    nosnoc.error('pds_rbp_not_one', "You are using an rk scheme with its right boundary point (c_n) not equal to one. Please choose another scheme e.g. RADAU_IIA.")
                end
                if opts.rk_representation == RKRepresentation.differential
                    nosnoc.error('pds_differential',"Differential representation without lifting is unsupported for gradient complementarity systems. Use integral or lifted differential representation.")
                end
                obj.dcs = nosnoc.dcs.PDSObjects(model);
                obj.dcs.generate_variables(opts);
                obj.dcs.generate_equations(opts);
                obj.discrete_time_problem = nosnoc.discrete_time_problem.PDSObjects(obj.dcs, opts);
                obj.discrete_time_problem.populate_problem();
              otherwise
                nosnoc.error('unknown_model', "Unknown model type.")
            end
        end

        function solve(obj, plugin)
            arguments
                obj
                plugin = 'reg_homotopy'
            end
            obj.discrete_time_problem.create_solver(obj.solver_opts, plugin);
            
            obj.stats = obj.discrete_time_problem.solve();
        end

        function ret = get(obj, field)
            opts = obj.opts;
            try
                warning off vdx:indexing:dot_reference_returns_vdx_var
                var = obj.discrete_time_problem.w.(field);
                warning on vdx:indexing:dot_reference_returns_vdx_var
            catch
                if strcmp(field, 'T_final')
                    nosnoc.warning("terminal_time_for_non_time_optimal_ocp",...
                        "You are trying to get the final time from a non-time-optimal problem. Instead returning the fixed time.")
                    ret = obj.discrete_time_problem.p.T.val;
                    return
                else
                    nosnoc.error('nonexistant_field', [char(field) ' is not a valid field for this OCP.']);
                end
                % TODO@anton print list of valid fields.
            end
            if var.depth == 3
                if opts.right_boundary_point_explicit
                    ret = var(:,:,obj.opts.n_s).res;
                else
                    ret = [var(0,0,obj.opts.n_s).res,...
                        var(1:opts.N_stages,1:opts.N_finite_elements(1),end).res];
                end
            else
                indexing(1:var.depth) = {':'};
                ret = var(indexing{:}).res;
            end
        end

        function set_param(obj, param, index, value)
            arguments
                obj
                param
                index cell = {}
                value double = 0;
            end
            if ~obj.discrete_time_problem.p.has_var(param);
                nosnoc.error('nonexistant_param', [char(param) ' does not exist as a parameter for this OCP.']);
            end
            warning off vdx:indexing:dot_reference_returns_vdx_var
            obj.discrete_time_problem.p.(param)(index{:}).val = value;
            warning on vdx:indexing:dot_reference_returns_vdx_var
        end

        function set_x0(obj, x0)
            obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).init = x0;
            obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).lb = x0;
            obj.discrete_time_problem.w.x(0,0,obj.opts.n_s).ub = x0;
        end

        function ret = get_full(obj, field)
            opts = obj.opts;
            try
                warning off vdx:indexing:dot_reference_returns_vdx_var
                var = obj.discrete_time_problem.w.(field);
                warning on vdx:indexing:dot_reference_returns_vdx_var
            catch
                nosnoc.error('nonexistant_field',[char(field) ' is not a valid field for this OCP.']);
                % TODO @anton print list of valid fields.
            end
            indexing(1:var.depth) = {':'};
            ret = var(indexing{:}).res;
        end

        function t_grid = get_time_grid(obj)
            opts = obj.opts;
            if opts.use_fesd
                h = obj.discrete_time_problem.w.h(:,:).res;
            else
                h = obj.discrete_time_problem.p.T().val/(sum(opts.N_finite_elements))*(ones(1, sum(opts.N_finite_elements)));
            end
            if opts.use_speed_of_time_variables
                if opts.local_speed_of_time_variable
                    sot = repelem(obj.get("sot"), opts.N_finite_elements);
                    h = sot.*h;
                else
                    sot = obj.get("sot");
                    h = sot*h;
                end
            end
            t_grid = cumsum([0, h]);
        end

        function t_grid_full = get_time_grid_full(obj)
            opts = obj.opts;
            if opts.use_fesd
                h = obj.discrete_time_problem.w.h(:,:).res;
            else
                h = obj.discrete_time_problem.p.T().val/(sum(opts.N_finite_elements))*(ones(1, sum(opts.N_finite_elements)));
            end
            t_grid_full = 0;
            if opts.use_speed_of_time_variables
                if opts.local_speed_of_time_variable
                    sot = repelem(obj.get("sot"), opts.N_finite_elements);
                    h = sot.*h;
                else
                    sot = obj.get("sot");
                    h = sot*h;
                end
            end
            for ii = 1:length(h)
                start = t_grid_full(end);
                for jj = 1:opts.n_s
                    t_grid_full = [t_grid_full; start + opts.c_rk(jj)*h(ii)];
                end
            end
        end

        function t_grid = get_control_grid(obj)
            t_grid = [0];
            for ii=1:obj.opts.N_stages
                if obj.opts.use_fesd
                    h_sum = sum(obj.discrete_time_problem.w.h(ii,:).res);
                else
                    h_sum = obj.discrete_time_problem.p.T().val/obj.opts.N_stages;
                end
                t_grid = [t_grid, t_grid(end)+h_sum];
            end
        end

        function objective = get_objective(obj)
            objective = obj.discrete_time_problem.f_result;
        end

        function set(obj, varname, field, indices, value)
            if ~obj.discrete_time_problem.w.has_var(varname)
                nosnoc.error('nonexistant_field', [char(varname) ' is not a valid field for this ocp.']);
            end
            var = obj.discrete_time_problem.w.(varname);
            var(indices{:}).(field) = value;
        end
        
        function generate_cpp_solver(obj, solver_dir, plugin)
            if ~exist('plugin', 'var')
                plugin = 'scholtes_ineq';
            end
            % add problem_options
            oldwd = pwd;
            solver_json = jsonencode(obj.opts, "PrettyPrint", true, "ConvertInfAndNaN", false);
            fid = fopen([solver_dir '/problem_options.json'], "w");
            fprintf(fid, solver_json);
            obj.discrete_time_problem.create_solver(obj.solver_opts, plugin);
            obj.discrete_time_problem.solver.generate_cpp_solver(solver_dir);
        end

        function set_initial_active_set(obj, active_set)
        % Set initial active set.
            arguments
                obj
                active_set(1,1) 
            end
            switch class(obj.model)
              case "nosnoc.model.Pss"
                if ~strcmp(class(active_set), "nosnoc.activeset.Pss")
                    nosnoc.error('type_mismatch', 'Wrong type of active set object passed');
                end
              case "nosnoc.model.Heaviside"
                error('not_implemented')
              case "nosnoc.model.Cls"
                error('not_implemented')
              case "nosnoc.model.Pds"
                if ~strcmp(class(active_set), "nosnoc.activeset.Pds")
                    nosnoc.error('type_mismatch', 'Wrong type of active set object passed');
                end
              case "nosnoc.model.PDSObjects"
                error('not_implemented')
              otherwise
                nosnoc.error('unknown_model', "Unknown model type.")
            end
            obj.active_set = active_set;
        end
        
        function do_shift_initialization(obj)
        % This method does a shift initialization by moving each control interval to the left by one.
        %
        % Warning:
        %    This is currently experimental and not guaranteed to work for all discretization settings.
            rbp = ~obj.opts.right_boundary_point_explicit;
            vars = obj.discrete_time_problem.w.get_vars();
            for var_idx=1:numel(vars)
                var = vars{var_idx};
                if var.depth == 0
                    continue
                end
                for ii=1:obj.opts.N_stages
                    next_ii = min(obj.opts.N_stages, ii+1);
                    if var.depth == 1
                        var(ii).init = var(next_ii).res;
                        continue
                    end
                    for jj=1:obj.opts.N_finite_elements(ii)
                        if var.depth == 2
                            var(ii,jj).init = var(next_ii,jj).res;
                            continue
                        end
                        for kk=1:(size(var.indices,3) - 1)
                            var(ii,jj,kk).init = var(next_ii,jj,kk).res;
                        end
                    end
                end
            end
        end

        function do_warmstart(obj)
        % This method warmstarts the next solve with the results of the previous solve.
            obj.discrete_time_problem.w.init = obj.discrete_time_problem.w.res;
        end
    end
end
 
