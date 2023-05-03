classdef NosnocSolver < handle
    properties
        model
        settings
        problem

        solver
        solver_initialization
    end

    methods
        function obj = NosnocSolver(model, settings)
            tic
            [obj.solver,obj.solver_initialization,obj.model,obj.settings] = create_nlp_nosnoc(model,settings);
            obj.problem = obj.model.problem;
            solver_generating_time = toc;
            if settings.print_level >=2
                fprintf('Solver generated in in %2.2f s. \n',solver_generating_time);
            end
        end

        function set(obj, type, val)
            if strcmp(type, 'w0')
                % check if val is size w0 then set it
            else
                ind = obj.problem.(strcat('ind_', type));
                flat_ind = sort([ind{:}]);

                if iscell(val)
                    % TODO do automatic fallback on smaller dimension cell arrays
                    if ndims(val) == 2 && size(val, 1) == obj.model.dims.N_stages && size(val,2) == 1
                        for ii=1:obj.model.dims.N_stages
                            for v=ind(ii,:,:)
                                if ~isempty(v) && length(v{1}) == length(val{ii})
                                    obj.solver_initialization.w0(v{1}) = val{ii}; 
                                end
                            end
                        end
                    elseif ndims(val) == 2 && size(val, 1) == obj.model.dimensions.N_stages && size(val, 2) == obj.model.dimeisons.N_fe
                    end
                else
                    if ndims(val) == 1
                        if length(val) == length(ind)
                            obj.solver_initialization.w0(flat_ind) = val;
                        end
                    else
                        error('nosnoc: set should be a cell array or a flat array')
                    end
                end
            end
        end
        
        function varargout = solve(obj)
            solver = obj.solver;
            model = obj.model;
            settings = obj.settings;
            solver_initialization = obj.solver_initialization
            [results,stats,solver_initialization] = homotopy_solver(solver,model,settings,solver_initialization);
            total_time = sum(stats.cpu_time);
            results = extract_results_from_solver(model,settings,results);

            % get 00 values
            % TODO: Make function for parameters 00
            lambda_00 = [];
            gamma_00 = [];
            p_vt_00 = [];
            n_vt_00  = [];
            gamma_d00 = [];
            delta_d00 = [];
            y_gap00 = [];
            if settings.dcs_mode == 'CLS'
                % TODO: reconsider this if 0th element has an impulse
                y_gap00 = model.f_c_fun(x0);
                if model.friction_exists
                    switch settings.friction_model
                        case 'Polyhedral'
                            v0 = x0(model.dimensions.n_q+1:end);
                            D_tangent_0 = model.D_tangent_fun(x0);
                            v_t0 = D_tangent_0'*v0;
                            for ii = 1:model.dimensions.n_contacts
                                ind_temp = model.dimensions.n_t*ii-(model.dimensions.n_t-1):model.dimensions.n_t*ii;
                                gamma_d00 = [gamma_d00;norm(v_t0(ind_temp))/model.dimensions.n_t];
                                delta_d00 = [delta_d00;D_tangent_0(:,ind_temp)'*v0+gamma_d00(ii)];
                            end
                        case 'Conic'
                            v0 = x0(model.dimensions.n_q+1:end);
                            v_t0 = model.J_tangent_fun(x0)'*v0;
                            for ii = 1:model.dimensions.n_contacts
                                ind_temp = model.dimensions.n_t*ii-(model.dimensions.n_t-1):model.dimensions.n_t*ii;
                                v_ti0 = v0(ind_temp);
                                gamma_00 = [gamma_00;norm(v_ti0)];
                                switch settings.conic_model_switch_handling
                                    case 'Plain'
                                        % no extra vars
                                    case {'Abs','Lp'}
                                        p_vt_00 = [p_vt_00;max(v_ti0,0)];
                                        n_vt_00 = [n_vt_00;max(-v_ti0,0)];
                                end
                            end
                    end
                end
            else
                lambda_00 = full(model.lambda00_fun(x0,model.p_global_val));
            end
            p_comp = [model.p_val;x0;lambda_00(:);y_gap00(:);gamma_00(:);gamma_d00(:);delta_d00(:);p_vt_00(:);n_vt_00(:)];
            complementarity_iter_ell_inf = full(comp_res(results.w_opt,p_comp));
            switch dcs_mode
                case 'Step'
                    temp = [results.alpha_opt_extended.*results.lambda_n_opt_extended,(1-results.alpha_opt_extended).*results.lambda_p_opt_extended];
                    complementarity_iter_ell_1 = sum(temp(:));
                case 'Stewart'
                    % TODO: considert cross comps as well in the inf norm
                    temp = [results.theta_opt_extended.*results.lam_opt_extended];
                    complementarity_iter_ell_1 = sum(temp(:));
                case 'CLS'
                    %         TODO?
            end

            stats.total_time  = total_time;
            fprintf('\n');
            fprintf('-----------------------------------------------------------------------------------------------\n');
            if settings.use_fesd
                fprintf( ['OCP with the FESD ' char(irk_scheme) ' in ' char(irk_representation) ' mode with %d RK-stages, %d finite elements and %d control intervals.\n'],n_s,N_finite_elements(1),N_stages);
            else
                fprintf( ['OCP with the Std ' char(irk_scheme) ' in ' char(irk_representation) ' mode with %d RK-stages, %d finite elements and %d control intervals.\n'],n_s,N_finite_elements(1),N_stages);
            end

            fprintf('---------------------------------------------- Stats summary--------------------------\n');
            if sum(stats.cpu_time) < 60
                fprintf('H. iters\t CPU Time (s)\t Max. CPU (s)/iter\tMin. CPU (s)/iter \tComp. res.\n');
                fprintf('%d\t\t\t\t%2.2f\t\t\t%2.2f\t\t\t\t%2.2f\t\t\t\t%2.2e\t\t\t\t%2.2e\n',stats.homotopy_iterations,sum(stats.cpu_time),max(stats.cpu_time),min(stats.cpu_time),complementarity_iter_ell_inf);
            else
                fprintf('H. iters\t CPU Time (m)\t Max. CPU (m)/iter\tMin. CPU (m)/iter \tComp. res.\n');
                fprintf('%d\t\t\t\t%2.2f\t\t%2.2f\t\t\t\t%2.2f\t\t\t\t\t%2.2e\t\t\t\t%2.2e \n',stats.homotopy_iterations,sum(stats.cpu_time)/60,max(stats.cpu_time)/60,min(stats.cpu_time)/60,complementarity_iter_ell_inf);
            end
            fprintf('\n--------------------------------------------------------------------------------------\n');
            if settings.time_optimal_problem
                T_opt = results.w_opt(model.ind_t_final);
                fprintf('Time optimal problem solved with T_opt: %2.4f.\n',T_opt);
                fprintf('\n--------------------------------------------------------------------------------------\n');
            end

            %% Output
            varargout{1} = results;
            varargout{2} = stats;
            varargout{3} = model;
            varargout{4} = settings;
            varargout{5} = solver_initialization;
        end 
    end
end

