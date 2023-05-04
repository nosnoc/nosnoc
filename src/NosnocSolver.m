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
                if length(obj.solver_initialization.w0) == length(val)
                    obj.solver_initialization.w0 = val;
                else
                    error("nosnoc: if initializing w0 all at once you need to provide a vector of corresponding size.")
                end
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

        
        function [results,stats] = solve(obj)
            import casadi.*
            solver = obj.solver;
            model = obj.model;
            settings = obj.settings;
            solver_initialization = obj.solver_initialization;
            
            comp_res = model.comp_res;
            nabla_J_fun = model.nabla_J_fun;

            % Initial conditions
            sigma_k = settings.sigma_0;
            w0 = solver_initialization.w0;
            p_val = obj.getInitialParameters(x0);

            % Initialize Stats struct
            stats = struct();
            stats.cpu_time = [];
            stats.cpu_time_total = 0;
            stats.sigma_k = sigma_k;
            stats.homotopy_iterations = [];
            stats.solver_stats = [];
            stats.objective = [];
            stats.complementarity_stats = [full(comp_res(w0, p_val))];

            
            % Initialize Results struct
            results = struct;
            results.W = [w0];

            lbw_h = solver_initialization.lbw; ubw_h = solver_initialization.ubw;
            lbw_h(model.ind_h) = model.h_k(1);
            ubw_h(model.ind_h) = model.h_k(1);

            % homotopy loop
            complementarity_iter = 1;
            ii = 0;

            if print_level >= 3
                %     fprintf('\niter\t\tsigma\t\tcompl_res\tobjective\tCPU time\tNLP iters\tstatus \t inf_pr \t inf_du \n')
                fprintf('\niter \t\t sigma \t\t compl_res \t\t objective \t\t inf_pr \t\t inf_du \t\t CPU time \t\t NLP iters \t\t status \n')
            end

            while (complementarity_iter) > comp_tol && ii < N_homotopy && (sigma_k > sigma_N || ii == 0)
                % homotopy parameter update
                if ii == 0
                    sigma_k = settings.sigma_0;
                else
                    if isequal(settings.homotopy_update_rule,'linear')
                        sigma_k = settings.homotopy_update_slope*sigma_k;
                    elseif isequal(settings.homotopy_update_rule,'superlinear')
                        sigma_k = max(settings.sigma_N,min(settings.homotopy_update_slope*sigma_k,sigma_k^settings.homotopy_update_exponent));
                    else
                        error('For the homotopy_update_rule please select ''linear'' or ''superlinear''.')
                    end
                end
                stats.sigma_k = [stats.sigma_k, sigma_k];
                p_val(1) = sigma_k;

                % Using multi-solver
                if iscell(solver)
                    tic
                    results = solver{ii+1}('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg,'p',p_val);
                    cpu_time_iter = toc ;
                    stats.solver_stats = [stats.solver_stats, solver{ii+1}.stats];
                else
                    tic
                    results = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg,'p',p_val);
                    cpu_time_iter = toc ;
                    stats.solver_stats = [stats.solver_stats, solver.stats];;
                end
                
                if isequal(stats.return_status,'Infeasible_Problem_Detected')
                    obj.printInfeasibility(results);
                end

                % update timing stats
                stats.cpu_time = [stats.cpu_time,cpu_time_iter];
                stats.cpu_time_total = stats.cpu_time_total + cpu_time_iter;

                % update results output.
                w_opt = full(results.x);
                results.W = [results.W,w_opt]; % all homotopy iterations
                
                w0 = w_opt;
                

                % update complementarity and objective stats
                complementarity_iter = full(comp_res(w_opt, p_val));
                stats.complementarity_stats = [stats.complementarity_stats;complementarity_iter];
                objective = full(model.problem.objective_fun(w_opt, p_val));
                stats.objective = [stats.objective, objective];
                
                % update counter
                ii = ii+1;
                % Verbose
                if print_level >= 3
                    obj.printNLPIterInfo(stats)
                end
            end

            % polish homotopy solution with fixed active set.
            % TODO fix this!
            if polishing_step
                [results] = polishing_homotopy_solution(model,settings,results,sigma_k);
                complementarity_iter = results.complementarity_iter;
                stats.complementarity_stats = [stats.complementarity_stats;complementarity_iter];
                W = [W,results.w_opt];
            end

            % number of iterations
            stats.homotopy_iterations = ii;

            results = extract_results_from_solver(model,settings,results);

            obj.printSolverStats(results,stats);
        end

        function printInfeasibility(obj, results)
            warning('nosnoc:homotopy_solver:NLP_infeasible', 'NLP infeasible: try different mpcc_mode or check problem functions.');
            if obj.settings.print_details_if_infeasible
                print_problem_details(results,obj.model,obj.solver_initialization,[]);
            end
            if obj.settings.pause_homotopy_solver_if_infeasible
                %             error('nosnoc: infeasible problem encounterd - stopping for debugging.')
                keyboard
            end
        end

        function p_val = getInitialParameters(obj)
            model = obj.model;
            settings = obj.settings;
            
            x0 = solver_initialization.w0(1:model.dimensions.n_x);
            lambda00 = [];
            gamma_00 = [];
            p_vt_00 = [];
            n_vt_00  = [];
            gamma_d00 = [];
            delta_d00 = [];
            y_gap00 = [];
            switch settings.dcs_mode
              case 'Stewart'
                lambda00 = full(model.lambda00_fun(x0, model.p_global_val));
              case 'Step'
                lambda00 = full(model.lambda00_fun(x0, model.p_global_val));
              case 'CLS'
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
            end
            p_val = [model.p_val(:);x0(:);lambda00(:);y_gap00(:);gamma_00(:);gamma_d00(:);delta_d00(:);p_vt_00(:);n_vt_00(:)];
        end

        function printNLPIterInfo(obj, stats)
            solver_stats = stats.solver_stats(end);
            ii = length(solver_stats);
            
            if strcmp(obj.settings.nlpsol, 'ipopt')
                inf_pr = solver_stats.iterations.inf_pr(end);
                inf_du = solver_stats.iterations.inf_du(end);
                fprintf('%d \t\t %6.2e \t\t %6.2e \t\t %6.3f \t\t %6.2e \t\t %6.2e \t\t %6.3f \t\t %d \t\t %s \n',...
                    ii, stats.sigma_k(end), stats.complementarity_stats(end), objective,inf_pr,inf_du, ...
                    cpu_time_iter, stats.iter_count, stats.return_status);
            elseif strcmp(settings.nlpsol, 'snopt')
                % TODO: Findout snopt prim du inf log!
                inf_pr = nan;
                inf_du = nan;
                fprintf('%d \t\t %6.2e \t\t %6.2e \t\t %6.3f \t\t %6.2e \t\t %6.2e \t\t %6.3f \t\t %d \t\t %s \n',...
                    ii, stats.sigma_k(end), complementarity_iter, objective,inf_pr,inf_du, ...
                    cpu_time_iter, solver_stats.secondary_return_status, solver_stats.return_status);
                warning('todo: add missing log information')
            end
        end

        function printSolverStats(obj, results, stats)
            model = obj.model;
            settings = obj.settings;
            
            comp_res = model.comp_res;


            fprintf('\n');
            fprintf('-----------------------------------------------------------------------------------------------\n');
            if settings.use_fesd
                fprintf( ['OCP with the FESD ' char(settings.irk_scheme) ' in ' char(settings.irk_representation) ' mode with %d RK-stages, %d finite elements and %d control intervals.\n'],...
                    n_s,N_finite_elements(1),N_stages);
            else
                fprintf( ['OCP with the Std ' char(settings.irk_scheme) ' in ' char(settings.irk_representation) ' mode with %d RK-stages, %d finite elements and %d control intervals.\n'],...
                    n_s,N_finite_elements(1),N_stages);
            end

            fprintf('---------------------------------------------- Stats summary--------------------------\n');
            if stats.cpu_time_total < 60
                fprintf('H. iters\t CPU Time (s)\t Max. CPU (s)/iter\tMin. CPU (s)/iter \tComp. res.\n');
                fprintf('%d\t\t\t\t%2.2f\t\t\t%2.2f\t\t\t\t%2.2f\t\t\t\t%2.2e\t\t\t\t%2.2e\n',...
                    stats.homotopy_iterations,stats.cpu_time_total,max(stats.cpu_time),min(stats.cpu_time),stats.complementarity_stats(end));
            else
                fprintf('H. iters\t CPU Time (m)\t Max. CPU (m)/iter\tMin. CPU (m)/iter \tComp. res.\n');
                fprintf('%d\t\t\t\t%2.2f\t\t%2.2f\t\t\t\t%2.2f\t\t\t\t\t%2.2e\t\t\t\t%2.2e \n',...
                    stats.homotopy_iterations,stats.cpu_time_total/60,max(stats.cpu_time)/60,min(stats.cpu_time)/60,stats.complementarity_stats(end));
            end
            fprintf('\n--------------------------------------------------------------------------------------\n');
            if settings.time_optimal_problem
                T_opt = results.w_opt(model.ind_t_final);
                fprintf('Time optimal problem solved with T_opt: %2.4f.\n',T_opt);
                fprintf('\n--------------------------------------------------------------------------------------\n');
            end
        end
        
    end
end

