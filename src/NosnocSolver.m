% BSD 2-Clause License

% Copyright (c) 2022, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% This file is part of NOSNOC.

classdef NosnocSolver < handle
    properties
        model     % Nosnoc Model
        settings  % Nosnoc Options
        problem   % Nosnoc reformulated NLP

        solver                 % CasADi solver function

        p_val % parameter values (homotopy parameters + problem parameters)
    end

    methods
        function obj = NosnocSolver(model, settings)
            import casadi.*
            tic
            model.verify_and_backfill(settings);
            model.generate_variables(settings);
            model.generate_equations(settings);

            % calculate homotopy
            if settings.N_homotopy == 0
                settings.N_homotopy = ceil(abs(log(settings.sigma_N / settings.sigma_0) / log(settings.homotopy_update_slope)));
                % TODO: compute
                if ~strcmp(settings.homotopy_update_rule, 'linear')
                    warning('computing N_homotopy automatically only supported for linear homotopy_update_rule');
                end
            end

            settings.preprocess();
            
            obj.model = model;
            obj.settings = settings;

            problem = NosnocProblem(settings, model.dims, model);
            obj.problem = problem;

            w = problem.w;
            g = problem.g;
            p = problem.p;
            if ~isempty(settings.ipopt_callback)
                settings.opts_casadi_nlp.iteration_callback = NosnocIpoptCallback('a_callback', model, problem, settings, length(w),length(g),length(p));
            end

            casadi_nlp = struct('f', problem.cost, 'x', w, 'g', g, 'p', p);

            % TODO: Possible issue raise to casadi: allow unknown fields in options passed
            if strcmp(settings.nlpsol, 'ipopt')
                opts_casadi_nlp = rmfield(settings.opts_casadi_nlp, 'snopt');
            elseif strcmp(settings.nlpsol, 'snopt')
                opts_casadi_nlp = rmfield(settings.opts_casadi_nlp, 'ipopt');
            end

            if ~settings.multiple_solvers
                solver = nlpsol(settings.solver_name, settings.nlpsol, casadi_nlp, opts_casadi_nlp);
            else
                solver = {};
                sigma_k = settings.sigma_0;
                for k = 1:settings.N_homotopy
                    opts_casadi_nlp.ipopt.mu_init = sigma_k * 1e-1;
                    opts_casadi_nlp.ipopt.mu_target = sigma_k * 1e-1;
                    opts_casadi_nlp.ipopt.bound_relax_factor = sigma_k^2 * 1e-2;
                    opts_casadi_nlp.ipopt.mu_strategy = 'monotone';
                    if k == 1
                        opts_casadi_nlp.ipopt.warm_start_init_point = 'yes';
                        opts_casadi_nlp.ipopt.warm_start_bound_push = 1e-4 * sigma_k;
                        opts_casadi_nlp.ipopt.warm_start_mult_bound_push = 1e-4 * sigma_k;
                    end
                    solver{k} = nlpsol(settings.solver_name, 'ipopt', casadi_nlp, opts_casadi_nlp);
                    % TODO: make homotopy update function and reuse here.
                    sigma_k = settings.homotopy_update_slope*sigma_k;
                end
            end
            obj.solver = solver;

            if ~isempty(settings.ipopt_callback)
                settings.opts_casadi_nlp.iteration_callback.solver = solver;
            end

            

            if settings.print_level > 5
                problem.print();
            end

            solver_generating_time = toc;
            if settings.print_level >=2
                fprintf('Solver generated in in %2.2f s. \n',solver_generating_time);
            end
        end

        function set(obj, type, val)
            if strcmp(type, 'w0')
                if length(obj.problem.w0) == length(val)
                    obj.problem.w0 = val;
                else
                    error("nosnoc: if initializing w0 all at once you need to provide a vector of corresponding size.")
                end
            else
                % This line uses the index sets collected during the creation of the NosnocProblem and automatically gets
                % the one of the form 'ind_<type>'. This makes this set generic for all variable types.
                ind = obj.problem.(strcat('ind_', type));
                if iscell(ind)
                    flat_ind = sort([ind{:}]);
                else
                    flat_ind = ind;
                end
                % TODO: Interpolation

                if iscell(val)
                    % If the passed value is an N_stage by 1 cell array we assume this initialization is done stage wise
                    if ismatrix(val) && size(val, 1) == obj.settings.N_stages && size(val,2) == 1
                        for ii=1:obj.settings.N_stages
                            % All variables of each stage are set to the same value
                            for v=ind(ii,:,:)
                                % NOTE: isempty check is needed for possibly unused rk-stage level cells (like in the case of rbp, etc.)
                                if ~isempty(v) && length(v{1}) == length(val{ii})
                                    obj.problem.w0(v{1}) = val{ii};
                                end
                            end
                        end
                    % Otherwise if we have an initialization of the form N_stages-by-N_fe we do the same but finite-element-wise
                    elseif ismatrix(val) && size(val, 1) == obj.settings.N_stages && size(val, 2) == obj.settings.N_finite_elements
                        for ii=1:obj.settings.N_stages
                            for jj=1:obj.settings.N_finite_elements
                                % All variables of each finite element are set to the same value
                                for v=ind(ii,jj,:)
                                    % NOTE: isempty check is needed for possibly unused rk-stage level cells (like in the case of rbp, cls, etc.)
                                    if ~isempty(v) && length(v{1}) == length(val{ii,jj})
                                        obj.problem.w0(v{1}) = val{ii,jj};
                                    end
                                end
                            end
                        end
                    else
                        error("nosnoc: Initialization either needs to be Stagewise or Finite Element wise, i.e. provide either an N_stage-by-1 or N_stage-by-N_fe cell array");
                    end
                    % Otherwise we assume that we are initializing via a flat array and we simply check for the same length
                else
                    if ismatrix(val) && length(val) == length(flat_ind)
                        obj.problem.w0(flat_ind) = val;
                    else
                        error('nosnoc: set should be a cell array or a flat array')
                    end
                end
            end
        end

        function [results,stats] = solve(obj)
            solver = obj.solver;
            model = obj.model;
            settings = obj.settings;
            problem = obj.problem;

            comp_res = problem.comp_res;

            % Initial conditions
            sigma_k = settings.sigma_0;
            w0 = problem.w0;
            lbw = problem.lbw; ubw = problem.ubw;
            lbg = problem.lbg; ubg = problem.ubg;
            obj.compute_initial_parameters();

            % Initialize Stats struct
            stats = struct();
            stats.cpu_time = [];
            stats.cpu_time_total = 0;
            stats.sigma_k = sigma_k;
            stats.homotopy_iterations = [];
            stats.solver_stats = [];
            stats.objective = [];
            stats.complementarity_stats = [full(comp_res(w0, obj.p_val))];


            % Initialize Results struct
            results = struct;
            results.W = w0;
            results.nlp_results = [];

            % homotopy loop
            complementarity_iter = 1;
            ii = 0;

            if settings.print_level >= 3
                %     fprintf('\niter\t\tsigma\t\tcompl_res\tobjective\tCPU time\tNLP iters\tstatus \t inf_pr \t inf_du \n')
                fprintf('\niter\t sigma \t\t compl_res\t inf_pr \t inf_du \t objective \t CPU time \t NLP iter\t status \n')
            end

            while (complementarity_iter) > settings.comp_tol && ii < settings.N_homotopy && (sigma_k > settings.sigma_N || ii == 0)
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
                obj.p_val(1) = sigma_k;

                % Using multi-solver
                if iscell(solver)
                    tic
                    nlp_results = solver{ii+1}('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg,'p',obj.p_val);
                    cpu_time_iter = toc ;
                    stats.solver_stats = [stats.solver_stats, solver{ii+1}.stats];
                else
                    tic
                    nlp_results = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg,'p',obj.p_val);
                    cpu_time_iter = toc ;
                    stats.solver_stats = [stats.solver_stats, solver.stats];
                end
                results.nlp_results = [results.nlp_results, nlp_results];

                if isequal(stats.solver_stats(end).return_status,'Infeasible_Problem_Detected')
                    obj.printInfeasibility(results);
                end

                % update timing stats
                stats.cpu_time = [stats.cpu_time,cpu_time_iter];
                stats.cpu_time_total = stats.cpu_time_total + cpu_time_iter;

                % update results output.
                w_opt = full(nlp_results.x);
                results.W = [results.W,w_opt]; % all homotopy iterations
                w0 = w_opt;

                % update complementarity and objective stats
                complementarity_iter = full(comp_res(w_opt, obj.p_val));
                stats.complementarity_stats = [stats.complementarity_stats;complementarity_iter];
                objective = full(obj.problem.objective_fun(w_opt, obj.p_val));
                stats.objective = [stats.objective, objective];

                % update counter
                ii = ii+1;
                % Verbose
                if settings.print_level >= 3
                    obj.printNLPIterInfo(stats)
                end
            end

            % polish homotopy solution with fixed active set.
            % TODO fix this!
            if settings.polishing_step
                [results] = polish_homotopy_solution(model,problem,settings,results,sigma_k);
                complementarity_iter = results.complementarity_iter;
                stats.complementarity_stats = [stats.complementarity_stats;complementarity_iter];
                W = [W,results.w];
            end

            % number of iterations
            stats.homotopy_iterations = ii;

            results = extract_results_from_solver(model,problem,settings,results);

            % check if solved to required accuracy
            stats.converged = obj.is_converged(stats);
            stats.constraint_violation = obj.compute_constraint_violation(results.w);

            obj.print_solver_stats(results,stats);
        end

        function violation = compute_constraint_violation(obj, w)
            problem = obj.problem;
            ubw_violation = max(max(w - problem.ubw), 0);
            lbw_violation = max(max(problem.lbw - w), 0);
            g_val = full(obj.problem.g_fun(w, obj.p_val));
            ubg_violation = max(max(g_val - problem.ubg), 0);
            lbg_violation = max(max(problem.lbg - g_val), 0);
            violation = max([lbg_violation, ubg_violation, lbw_violation, ubw_violation]);
        end

        function converged = is_converged(obj, stats)
            if strcmp(obj.settings.nlpsol, 'ipopt')
                last_stats = stats.solver_stats(end);
                converged = 0;
                if isfield(last_stats, 'iterations')
                    inf_pr = last_stats.iterations.inf_pr(end);
                    inf_du = last_stats.iterations.inf_du(end);
                    if inf_pr < obj.settings.opts_casadi_nlp.ipopt.tol && inf_du < obj.settings.opts_casadi_nlp.ipopt.tol ...
                            && stats.complementarity_stats(end) < 10 * obj.settings.comp_tol
                        converged = 1;
                    end
                else
                    converged = 0;
                end
            else
                % TODO..
                converged = [];
            end
        end


        function printInfeasibility(obj, results)
            % warning('nosnoc:homotopy_solver:NLP_infeasible', 'NLP infeasible: try different mpcc_mode or check problem functions.');
            if obj.settings.print_details_if_infeasible
                print_problem_details(results,obj.model,obj.problem, []);
            end
            if obj.settings.pause_homotopy_solver_if_infeasible
                %             error('nosnoc: infeasible problem encounterd - stopping for debugging.')
                keyboard
            end
        end


        function compute_initial_parameters(obj)
            model = obj.model;
            settings = obj.settings;

            x0 = obj.problem.w0(1:model.dims.n_x);
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
                y_gap00 = max(0, model.f_c_fun(x0));
                if model.friction_exists
                    switch settings.friction_model
                      case 'Polyhedral'
                        v0 = x0(model.dims.n_q+1:end);
                        D_tangent_0 = model.D_tangent_fun(x0);
                        v_t0 = D_tangent_0'*v0;
                        for ii = 1:model.dims.n_contacts
                            ind_temp = model.dims.n_t*ii-(model.dims.n_t-1):model.dims.n_t*ii;
                            gamma_d00 = [gamma_d00;norm(v_t0(ind_temp))/model.dims.n_t];
                            delta_d00 = [delta_d00;D_tangent_0(:,ind_temp)'*v0+gamma_d00(ii)];
                        end
                      case 'Conic'
                        v0 = x0(model.dims.n_q+1:end);
                        v_t0 = model.J_tangent_fun(x0)'*v0;
                        for ii = 1:model.dims.n_contacts
                            ind_temp = model.dims.n_t*ii-(model.dims.n_t-1):model.dims.n_t*ii;
                            v_ti0 = v0(ind_temp);
                            gamma_00 = [gamma_00;norm(v_ti0)];
                            switch settings.conic_model_switch_handling
                              case 'Plain'
                                % no extra vars
                              case {'Abs','Lp'}
                                p_vt_00 = [p_vt_00; max(v_ti0,0)];
                                n_vt_00 = [n_vt_00; max(-v_ti0,0)];
                            end
                        end
                    end
                end
            end
            obj.p_val = [obj.problem.p0(:);x0(:);lambda00(:);y_gap00(:);gamma_00(:);gamma_d00(:);delta_d00(:);p_vt_00(:);n_vt_00(:)];
        end

        function printNLPIterInfo(obj, stats)
            solver_stats = stats.solver_stats(end);
            ii = size(stats.solver_stats, 2);

            if strcmp(obj.settings.nlpsol, 'ipopt')
                if isfield(solver_stats, 'iterations')
                    inf_pr = solver_stats.iterations.inf_pr(end);
                    inf_du = solver_stats.iterations.inf_du(end);
                    objective = solver_stats.iterations.obj(end);
                else
                    inf_pr = nan;
                    inf_du = nan;
                    objective = nan;
                end
                fprintf('%d\t%6.2e\t %6.2e\t %6.2e\t %6.2e \t %6.2e \t %6.3f \t %d \t %s \n',...
                    ii, stats.sigma_k(end), stats.complementarity_stats(end), inf_pr,inf_du, ...
                    objective, stats.cpu_time(end), solver_stats.iter_count, solver_stats.return_status);
            elseif strcmp(obj.settings.nlpsol, 'snopt')
                % TODO: Findout snopt prim du inf log!
                inf_pr = nan;
                inf_du = nan;
                fprintf('%d\t%6.2e\t %6.2e\t %6.2e \t %6.2e \t %6.2e \t %6.3f \t %s \t %s \n',...
                    ii, stats.sigma_k(end), stats.complementarity_stats(end), inf_pr,inf_du, ...
                    stats.objective(end), stats.cpu_time(end), solver_stats.secondary_return_status, solver_stats.return_status);
                % TODO: add missing log information
            end
        end

        function print_iterate(obj, iterate, only_violations, filename)
            if exist('filename')
                delete(filename);
                fileID = fopen(filename, 'w');
            else
                fileID = 1;
            end
            if ~exist('only_violations')
                only_violations = 0;
            end
            fileID = 1;
            fprintf(fileID, "\nw\t\t\tlbw\t\tubw\titerate\n");
            for i = 1:length(obj.problem.lbw)
                if ~only_violations || (iterate(i) < obj.problem.lbw(i) || iterate(i) > obj.problem.ubw(i))
                    expr_str = pad(formattedDisplayText(obj.problem.w(i)), 20);
                    lb_str = pad(sprintf('%.2e', obj.problem.lbw(i)), 10);
                    ub_str = pad(sprintf('%.2e', obj.problem.ubw(i)), 10);
                    iterate_str = pad(sprintf('%.2e', iterate(i)), 10);
                    fprintf(fileID, "%s\t%s\t%s\t%s\n", expr_str, lb_str, ub_str, iterate_str);
                end
            end

            % constraints
            g_val = full(obj.problem.g_fun(iterate, obj.p_val));
            fprintf(fileID, "\ni\tlbg\t\t ubg\t\t g_val\t\tg_expr\n");
            for i = 1:length(obj.problem.lbg)
                if ~only_violations || (g_val(i) < obj.problem.lbg(i) || g_val(i) > obj.problem.ubg(i))
                    expr_str = formattedDisplayText(obj.problem.g(i));
                    lb_str = pad(sprintf('%.2e', obj.problem.lbg(i)), 12);
                    ub_str = pad(sprintf('%.2e', obj.problem.ubg(i)), 12);
                    fprintf(fileID, "%d\t%s\t%s\t%.2e\t%s\n", i, lb_str, ub_str, g_val(i), expr_str);
                end
            end

        end
        function print_solver_stats(obj, results, stats)
            % model = obj.model;
            % dims = model.dims;
            % settings = obj.settings;

            % fprintf('\n---------------------------------------------- Stats summary--------------------------\n');
            % if stats.cpu_time_total < 60
            %     fprintf('H. iters\t CPU Time (s)\t Max. CPU (s)/iter\tMin. CPU (s)/iter \tComp. res.\n');
            %     fprintf('%d\t\t\t\t%2.2f\t\t\t%2.2f\t\t\t\t%2.2f\t\t\t\t%2.2e\t\t\t\t%2.2e\n',...
            %         stats.homotopy_iterations, stats.cpu_time_total, max(stats.cpu_time),min(stats.cpu_time), stats.complementarity_stats(end));
            % else
            %     fprintf('H. iters\t CPU Time (m)\t Max. CPU (m)/iter\tMin. CPU (m)/iter \tComp. res.\n');
            %     fprintf('%d\t\t\t\t%2.2f\t\t%2.2f\t\t\t\t%2.2f\t\t\t\t\t%2.2e\t\t\t\t%2.2e \n',...
            %         stats.homotopy_iterations,stats.cpu_time_total/60, max(stats.cpu_time)/60, min(stats.cpu_time)/60, stats.complementarity_stats(end));
            % end
            % fprintf('\n');
        end

    end
end
