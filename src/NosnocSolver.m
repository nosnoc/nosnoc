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
        solver_initialization  % initialization for the generated NLP
    end

    methods
        function obj = NosnocSolver(model, settings)
            tic
            import casadi.*;

            [model,settings] = model_reformulation_nosnoc(model,settings); % TODO Move this outside solver to NosnocModel Class

            settings.create_butcher_tableu(model); % TODO this should live somewhere else. (i.e. butcher tableu should not be in settings)

            obj.model = model;
            obj.settings = settings;

            problem = NosnocProblem(settings, model.dims, model);
            obj.problem = problem;

            w = problem.w;
            g = problem.g;
            p = problem.p;
            J_fun = problem.cost_fun;
            comp_res = problem.comp_res;
            comp_res_fesd = problem.comp_fesd;
            comp_res_std = problem.comp_std;

            casadi_nlp = struct('f', problem.cost, 'x', w, 'g', g,'p',p);

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
                    opts_casadi_nlp = settings.opts_casadi_nlp;
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

            %% Define CasADi function for the switch indicator function.
            nu_fun = Function('nu_fun', {w,p},{problem.nu_vector});

            % TODO Clean this up
            % I.e: - No longer duplicate things in model and solver.
            %      - Separate model generation.
            obj.model.problem = problem;
            obj.model.g = g;
            obj.model.w = w;
            obj.model.p = p;
            obj.model.J = problem.cost;
            obj.model.J_fun = J_fun;
            obj.model.comp_res = comp_res;
            obj.model.comp_res_fesd = comp_res_fesd;
            obj.model.comp_res_std = comp_res_std;
            obj.model.nu_fun = nu_fun;

            % create CasADi function for objective gradient.
            nabla_J = problem.cost.jacobian(obj.model.w);
            nabla_J_fun = Function('nabla_J_fun', {w,p},{nabla_J});
            obj.model.nabla_J = nabla_J;
            obj.model.nabla_J_fun = nabla_J_fun;

            if settings.print_level > 5
                problem.print();
            end

            %% Model update: all index sets and dims
            % TODO: Maybe just return the problem, currently trying not to break compatibility for now.
            obj.model.ind_x = [problem.ind_x0.'; flatten_ind(problem.ind_x)];
            obj.model.ind_v = sort(flatten_ind(problem.ind_v));
            obj.model.ind_z_all = problem.ind_z_all; %TODO fix this by breaking compat
            obj.model.ind_u = problem.ind_u;
            obj.model.ind_h = flatten_ind(problem.ind_h);
            obj.model.ind_sot = flatten_ind(problem.ind_sot);
            obj.model.ind_t_final  = problem.ind_t_final;
            obj.model.p_val = problem.p0;

            %% Store solver initialization data
            solver_initialization.w0 = problem.w0;
            solver_initialization.lbw = problem.lbw;
            solver_initialization.ubw = problem.ubw;
            solver_initialization.lbg = problem.lbg;
            solver_initialization.ubg = problem.ubg;

            obj.solver_initialization = solver_initialization;
            solver_generating_time = toc;
            if settings.print_level >=2
                fprintf('Solver generated in in %2.2f s. \n',solver_generating_time);
            end
        end

        function set(obj, type, val)
            if strcmp(type, 'w0')
                if length(obj.solver_initialization.w0) == length(val)
                    obj.solver_initialization.w0 = val;
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

                if iscell(val)
                    % If the passed value is an N_stage by 1 cell array we assume this initialization is done stage wise
                    if ndims(val) == 2 && size(val, 1) == obj.model.dims.N_stages && size(val,2) == 1
                        for ii=1:obj.model.dims.N_stages
                            % All variables of each stage are set to the same value
                            % TODO: Interpolation
                            for v=ind(ii,:,:)
                                % NOTE: isempty check is needed for possibly unused rk-stage level cells (like in the case of rbp, etc.)
                                if ~isempty(v) && length(v{1}) == length(val{ii})
                                    obj.solver_initialization.w0(v{1}) = val{ii};
                                end
                            end
                        end
                    % Otherwise if we have an initialization of the form N_stages-by-N_fe we do the same but finite-element-wise
                    elseif ndims(val) == 2 && size(val, 1) == obj.model.dims.N_stages && size(val, 2) == obj.model.dimeisons.N_fe
                        for ii=1:obj.model.dims.N_stages
                            for jj=1:obj.model.dims.N_fe
                                % All variables of each finite element are set to the same value
                                % TODO: Interpolation
                                for v=ind(ii,jj,:)
                                    % NOTE: isempty check is needed for possibly unused rk-stage level cells (like in the case of rbp, cls, etc.)
                                    if ~isempty(v) && length(v{1}) == length(val{ii,jj})
                                        obj.solver_initialization.w0(v{1}) = val{ii,jj};
                                    end
                                end
                            end
                        end
                    else
                        error("nosnoc: Initialization either needs to be Stagewise or Finite Element wise, i.e. provide either an N_stage-by-1 or N_stage-by-N_fe cell array");
                    end
                    % Otherwise we assume that we are initializing via a flat array and we simply check for the same length
                else
                    if ndims(val) == 2 && length(val) == length(ind)
                        obj.solver_initialization.w0(flat_ind) = val;
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
            lbw = solver_initialization.lbw; ubw = solver_initialization.ubw;
            lbg = solver_initialization.lbg; ubg = solver_initialization.ubg;
            p_val = obj.getInitialParameters();

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
            results.nlp_results = [];

            % homotopy loop
            complementarity_iter = 1;
            ii = 0;

            if settings.print_level >= 3
                %     fprintf('\niter\t\tsigma\t\tcompl_res\tobjective\tCPU time\tNLP iters\tstatus \t inf_pr \t inf_du \n')
                fprintf('\niter \t\t sigma \t\t compl_res \t\t objective \t\t inf_pr \t\t inf_du \t\t CPU time \t\t NLP iters \t\t status \n')
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
                p_val(1) = sigma_k;

                % Using multi-solver
                if iscell(solver)
                    tic
                    nlp_results = solver{ii+1}('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg,'p',p_val);
                    cpu_time_iter = toc ;
                    stats.solver_stats = [stats.solver_stats, solver{ii+1}.stats];
                else
                    tic
                    nlp_results = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg,'p',p_val);
                    cpu_time_iter = toc ;
                    stats.solver_stats = [stats.solver_stats, solver.stats];;
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
                complementarity_iter = full(comp_res(w_opt, p_val));
                stats.complementarity_stats = [stats.complementarity_stats;complementarity_iter];
                objective = full(model.problem.objective_fun(w_opt, p_val));
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
                [results] = polish_homotopy_solution(model,settings,results,sigma_k);
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

            x0 = obj.solver_initialization.w0(1:model.dims.n_x);
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
            ii = size(solver_stats, 2);

            if strcmp(obj.settings.nlpsol, 'ipopt')
                if isfield(solver_stats, 'iterations')
                    inf_pr = solver_stats.iterations.inf_pr(end);
                    inf_du = solver_stats.iterations.inf_du(end);
                else
                    inf_pr = nan;
                    inf_du = nan;
                end
                
                fprintf('%d \t\t %6.2e \t\t %6.2e \t\t %6.3f \t\t %6.2e \t\t %6.2e \t\t %6.3f \t\t %d \t\t %s \n',...
                    ii, stats.sigma_k(end), stats.complementarity_stats(end), stats.objective(end),inf_pr,inf_du, ...
                    stats.cpu_time(end), solver_stats.iter_count, solver_stats.return_status);
            elseif strcmp(obj.settings.nlpsol, 'snopt')
                % TODO: Findout snopt prim du inf log!
                inf_pr = nan;
                inf_du = nan;
                fprintf('%d \t\t %6.2e \t\t %6.2e \t\t %6.3f \t\t %6.2e \t\t %6.2e \t\t %6.3f \t\t %s \t\t %s \n',...
                    ii, stats.sigma_k(end), stats.complementarity_stats(end), stats.objective(end),inf_pr,inf_du, ...
                    stats.cpu_time(end), solver_stats.secondary_return_status, solver_stats.return_status);
                % TODO: add missing log information
            end
        end

        function printSolverStats(obj, results, stats)
            model = obj.model;
            dims = model.dims;
            settings = obj.settings;

            comp_res = model.comp_res;


            fprintf('\n');
            fprintf('-----------------------------------------------------------------------------------------------\n');
            if settings.use_fesd
                fprintf( ['OCP with the FESD ' char(settings.irk_scheme) ' in ' char(settings.irk_representation) ' mode with %d RK-stages, %d finite elements and %d control intervals.\n'],...
                    dims.n_s,dims.N_finite_elements(1),dims.N_stages);
            else
                fprintf( ['OCP with the Std ' char(settings.irk_scheme) ' in ' char(settings.irk_representation) ' mode with %d RK-stages, %d finite elements and %d control intervals.\n'],...
                    dims.n_s,dims.N_finite_elements(1),dims.N_stages);
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
