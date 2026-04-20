classdef Rtmpc < nosnoc.mpc.Base
    % Real-time MPC based on the sequential quadratic programing with complementarity constraints.
    % Inputs for constructor:
    %  - model (required): system model object with dynamics, constraints, and cost definition.
    %  - mpc_options, problem_options, solver_options, qpec_options, qp_options (optional): structs configuring MPC setup, problem formulation, solver, QPCC handling, and QP solver settings.


    properties %(Access=private)
        mpc_options nosnoc.mpc.Options % nosnoc mpc specific options.

        % nosnoc ocp solver stuff (for stroing all data, but creatse also an mpcc solver - may not be needed at al ater point)
        model % nosnoc model.
        problem_options nosnoc.Options % nosnoc optimal control problem options.
        solver_options % nosnoc.reg_homotopy.Options % options for ocp_solver, this creates a nosnoc solver where all data is stored (linearization points, parameters, functions); (may be later split)
        ocp_solver nosnoc.ocp.Solver % nosnoc OCP solver used to generate feedback

        % qpec stuff
        qpec_functions % all qpec functions (residuals, first and second order derivatives) (TODO: should live in qpec class later)
        qpec % Qpec object
        qpec_options % qpec solver options: Can be homotopy_options if mpccsol used, mpecopt options if mpecopt, lcqpow optiosn, gurobi options
        % qp probing/feedback
        qp_options
        qp_solver
        qp_acceptable_status
        % book-keeping
        % last_solve_successful(1,1) logical = false % true if last solve was successful, false otherwise. % TODO REMPVE tHIS
        feedback_phase_successful(1,1) logical = true % true if feedback QPCC/QP was solved; otherwise use some fallbakc strategy
        very_first_mpc_step(1,1) logical = true % flag if the mpc was just started. used to solve first problem fully (if this option is turned on)
        stats
    end

    methods

        function obj = Rtmpc(model, mpc_options, problem_options, solver_options, qpec_options, qp_options)
            % Constructor: initialize SQPCC-MPC object, create OCP solver and QPCC structures.

            obj.mpc_options = mpc_options; % MPC-level options (shift, warmstart, etc.)
            % ----- handle optional arguments -----
            if nargin < 4 || isempty(solver_options)
                % default MPCC options if none provided
                solver_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosonc problems.
                solver_options.complementarity_tol = 1e-8; % Value to drive the complementarity residual to.
                solver_options.N_homotopy = 8; % Maximum number of homotopy iterations.
                solver_options.print_level = 0;
                solver_options.sigma_0 = 1e0;
                solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; %
                solver_options.homotopy_steering_strategy = "DIRECT";
                solver_options.lift_complementarities = 0;
            end

            if nargin < 5 || isempty(qpec_options)
                % Default qpec solver options, reg homotopy;
                qpec_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosonc problems.
                qpec_options.complementarity_tol = 1e-8; % Value to drive the complementarity residual to.
                qpec_options.N_homotopy = 8; % Maximum number of homotopy iterations.
                qpec_options.print_level = 0;
                qpec_options.sigma_0 = 1e0;
                qpec_options.homotopy_steering_strategy = "DIRECT";
                qpec_options.lift_complementarities = 0;

                % Nlp solver hints
                qpec_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; %
                qpec_options.opts_casadi_nlp.ipopt.hessian_constant = 'yes';
                qpec_options.opts_casadi_nlp.ipopt.jac_c_constant   = 'yes';
                qpec_options.opts_casadi_nlp.ipopt.max_iter = 200;
                qpec_options.opts_casadi_nlp.ipopt.nlp_scaling_method = 'none';
                % qpec_options.opts_casadi_nlp.ipopt.mu_strategy = 'monotone';
                % qpec_options.opts_casadi_nlp.ipopt.mehrotra_algorithm = 'yes';
                qpec_options.opts_casadi_nlp.ipopt.accept_every_trial_step = 'yes';
                qpec_options.opts_casadi_nlp.ipopt.mu_init = 1e-1;
            end
            if nargin < 6 || isempty(qp_options)
                % if qp probing is active, default option is ipopt
                if mpc_options.use_probing_qp || mpc_options.use_feedback_qp
                    % TODO: consider putting this contruction inside the QPCC (should we pass qp options to the rtmpc\qpce
                    % constructors?)
                    % create qp solver directly via the gurobi solver class (can be made also via qpec like ipopt below)
                    % obj.qp_options = nosnoc.qpec.GurobiOptions();
                    % obj.qp_options.method = 'qp';
                    % obj.qp_solver = nosnoc.qpec.GurobiQpecSolver(obj.qp_options); % This solver expect a qpec as input
                    % % create a QP, a Qpec class object, and build an ipopt qp
                    % % solver
                    % Default QP solver options (TODO; make an options file?)
                    obj.qp_options = struct;%  ipopt for a qp solver

                    obj.qp_options.ipopt.tol = 1e-8;
                    obj.qp_options.ipopt.max_iter = 500;
                    obj.qp_options.ipopt.tiny_step_tol = 1e-20;
                    obj.qp_options.ipopt.fixed_variable_treatment = 'make_constraint';
                    % QP structure hints (cant be used if QP is fully parametric?)
                    % obj.qp_options.ipopt.hessian_constant = 'yes';
                    % obj.qp_options.ipopt.jac_c_constant   = 'yes';
                    % obj.qp_options.ipopt.jac_d_constant   = 'yes';
                    % obj.qp_options.ipopt.grad_f_constant  = 'yes';
                    % compling
                    % obj.qp_options.jit = true;
                    % obj.qp_options.compiler = 'shell';
                    % obj.qp_options.jit_options.compiler = 'gcc';
                    % obj.qp_options.jit_options.flags = {'-O2'};
                    % Mehrotra predictor–corrector and related
                    obj.qp_options.ipopt.mehrotra_algorithm = 'yes';
                    obj.qp_options.ipopt.accept_every_trial_step = 'yes';
                    obj.qp_options.ipopt.mu_init = 1e0;
                    % Solver and printing
                    obj.qp_options.ipopt.linear_solver = 'mumps';
                    obj.qp_options.ipopt.print_level = 0;
                    obj.qp_options.ipopt.print_timing_statistics = 'no';
                    % qp = obj.setup_qp_solver('qp_ipopt'); % Create QP object (storse data, sparsity and solver function, the qp isnt a propery of rtmpc class as necessary info are stored in the qpec object)
                    % obj.qp_solver = qp.solver;

                    obj.qp_acceptable_status = {'Solve_Succeeded', ...
                        'Solved_To_Acceptable_Level', ...
                        'Search_Direction_Becomes_Too_Small', ...
                        'Maximum_Iterations_Exceeded'};
                else
                    obj.qp_options = []; % no probing QP
                    obj.qp_solver = [];   % not used
                    obj.qp_acceptable_status = {};
                end
            end

            if ~obj.mpc_options.solve_advanced_problem
                obj.mpc_options.use_feedback_qp = false; % only qp without advanded step makes little sense; try `use_probing_qp` instead;
            end

            % All OCP Data
            obj.problem_options = problem_options; % Problem definition options (N, constraints, bounds)
            obj.solver_options = solver_options;   % Solver configuration for OCP (e.g. tolerances, NLP solver)
            obj.ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options); % Initialize and solve base OCP once

            % QPCC/QP specific
            obj.qpec_options = qpec_options; % QPCC solver options (for mpccsol, mpecopt, lcqpow, or gurobi)
            obj.qpec_functions = obj.build_qpec_symbolic_functions(); % Generate CasADi functions for QPCC (grad/Hess/constraints)
            obj.qpec = obj.setup_qpec_solver();               % Create QPCC object and qpec solver (storse data, sparsity and solver function)
            obj.setup_qp_solver();
            obj.initialize_qpec();                            % Initialize QPCC object and update bounds, and warm start

            % obj.last_solve_successful = true;   % since the constructor OCP solve succeeded (todo read it from ocp solveR)
            obj.very_first_mpc_step = true;     % only used to skip one re-solve in feedback

            obj.stats = struct(); % Initialize stats structure for logging feedback and preparation phases

            % Timing stats;
            obj.stats.feedback_times = {};
            obj.stats.preparation_times = {};
            % TODO: Add more detailed timings;

        end

        % ----------------------------------- MPC functions preparation and feedback ---------------------------

        function [stats] = do_preparation(obj,x0_pred)
            % This method should be called after get feedback to do the work of setting up the warm starting.
            % It does all the derivative evalautions to set up a new QPCC.

            preparation_timer = tic;
            preparation_solver_time = 0.0;
            if obj.mpc_options.do_shift_initialization
                obj.ocp_solver.do_shift_initialization();
                % Note: this shifts w0, and so it does implicitly y0, when we eval G and H to compute the y0
            else
                obj.ocp_solver.do_warmstart();
            end

            if obj.mpc_options.solve_advanced_problem
                fprintf('nosnoc mpc: Solving advanced problem in preparation phase. \n')
                if ~exist('x0_pred')
                    if obj.feedback_phase_successful
                        x0_pred = obj.get_predicted_state(); % shift from OCP solution;
                    else
                        x0_pred = obj.get_current_state(); % if the previous feedback phase failed; the prediction is probbably bad;
                    end
                end

                switch obj.mpc_options.advanced_problem_type
                    case 'full'
                        fprintf('nosnoc mpc: Solving advanced MPCC to convergence in preparation phase. \n')
                        obj.ocp_solver.set_x0(x0_pred); % TODO: implement x0_pred to be passed from outside, e.g. from simulatuon
                        obj.ocp_solver.solve(); % Solve advanced problem via mpccsol to convergence!
                        if isa(obj.solver_options,'nosnoc.ccopt.Options')
                            preparation_solver_time = obj.ocp_solver.stats.ccopt.total_wall_time;
                        else
                            preparation_solver_time = obj.ocp_solver.stats.wall_time_total;
                        end

                    case 'sqpec'
                        fprintf('nosnoc mpc: Solving advanced QPCC in preparation phase. \n')
                        obj.ocp_solver.set_x0(x0_pred);
                        obj.update_initial_value_in_qpec(); % updating in qpec, assumes it is updated in the ocp_solver, as it stores all linearizaiton points;
                        for ii = 1:obj.mpc_options.advanced_n_qpecs
                            stats.update_stats = obj.update_qpec_data();  % Evaluate and update QPCC matrices
                            results_qpec = obj.qpec.solve();
                            % stats.qpec_solved = 1;
                            p_val = obj.ocp_solver.discrete_time_problem.p.val; % problem parameters
                            w_k = obj.ocp_solver.discrete_time_problem.w.res;   % read linearization point from ocp solver

                            d_qpec    = full(results_qpec.x);
                            lam_g_new = full(results_qpec.lam_g);
                            lam_x_new = full(results_qpec.lam_x);
                            w_k       = w_k + d_qpec;

                            if isa(obj.qpec_options,'nosnoc.ccopt.Options')
                                preparation_solver_time = preparation_solver_time + obj.qpec.stats.ccopt.total_wall_time;
                            else
                                preparation_solver_time = preparation_solver_time + obj.qpec.stats.wall_time_total;
                            end
                        end
                        % Update values in OCP solver memory
                        obj.ocp_solver.discrete_time_problem.w.res = w_k;
                        obj.ocp_solver.discrete_time_problem.w.mult = lam_x_new;
                        obj.ocp_solver.discrete_time_problem.g.mult = lam_g_new;
                    otherwise
                        % none yet;
                end
            end
            % Prepare the QPCC for the feedback phase:
            fprintf('nosnoc mpc: Preparing QPCC in feedback phase. \n')
            % TODO: this update can be optional in advanced step with qpecs (effectivly removes corrector part)
            stats.update_stats = obj.update_qpec_data();  % Evaluate and update QPCC matrices
            % end
            stats.preparation_time = toc(preparation_timer);
            if obj.mpc_options.solve_advanced_problem
                stats.preparation_time = preparation_solver_time;
            end

            if isfield(obj.stats, 'preparation_stats')
                obj.stats.preparation_stats(end+1) = stats;
            else
                obj.stats.preparation_stats = stats;
            end
            obj.stats.preparation_times{end+1} = stats.preparation_time;
        end

        function [u, stats] = get_feedback(obj, x0)
            % Compute MPC control action using SQPCC step or full OCP solve if needed.
            feedback_timer = tic;
            stats.qp_solved = 0;
            stats.qpec_solved = 0;
            stats.qp_accepted = 0;

            if obj.very_first_mpc_step
                % --- first feedback: use already solved OCP, no QPCC yet ---
                fprintf('nosnoc mpc: First MPC iteration, using initial OCP solution.\n');
                u_res = obj.ocp_solver.get('u');
                u = u_res(:,1);
                stats.d_qpec_norm = 0;
                stats.feedback_time = toc(feedback_timer);
                obj.very_first_mpc_step = false;
                return
            end

            % Read latest solution (lin. point) and parameters
            p_val = obj.ocp_solver.discrete_time_problem.p.val; % problem parameters
            w_k = obj.ocp_solver.discrete_time_problem.w.res;   % read linearization point from ocp solver
            % TODO: obj.qpec_initialization.p(end-length(p_val)+1:end) = p_val;  % Update parameters of OCP if they are changed in the feedback phase (e.g. reference is updated)

            % Update inital value in ocp solver and in QPCC/QP.
            obj.ocp_solver.set_x0(x0);
            obj.update_initial_value_in_qpec(); % updating in qpec, assumes it is updated in the ocp_solver, as it stores all linearizaiton points;

            % some homotopy warm start option for reg homotopy (e.g. start with smaller homotopy parameter)
            if isa(obj.qpec_options,'nosnoc.reg_homotopy.Options')
                if obj.mpc_options.warmstart_qpec
                    if ~obj.very_first_mpc_step
                        obj.qpec.solver_opts.sigma_0 = obj.mpc_options.fast_sigma_0;
                        obj.qpec.w0 = obj.ocp_solver.stats.homotopy_iterations(:,1);
                    end
                end
            end

            if ~(obj.mpc_options.use_probing_qp || obj.mpc_options.use_feedback_qp)
                % Solve QPCC
                fprintf('nosnoc mpc: Solving QPCC in feedback phase. \n')
                results_qpec = obj.qpec.solve();
                stats.qpec_solved = obj.qpec.stats.success;
                stats.qp_solved = 0;
                f_k_new = full(obj.qpec_functions.f_fun(w_k+full(results_qpec.x), p_val));
                fprintf('QPCC objective %2.2f \n',f_k_new);
            else
                fprintf('nosnoc mpc: Solving QP in feedback phase. \n')
                [results_qp, stats_qp] = obj.qpec.solve_qp(); % via qpec class
                f_k_new = full(obj.qpec_functions.f_fun(w_k+full(results_qp.x), p_val));
                fprintf('QP objective %2.2f \n',f_k_new);
                % Check is QP solve succesful:
                if stats_qp.success
                    stats.qpec_solved = 0;
                    stats.qp_solved = 1;
                    % if ismember(stats_qp.return_status, obj.qp_acceptable_status)
                    %     stats.qp_solved =  1;
                    % end
                end
                % directly gurobi solver class interface
                % [results_qp, stats_qp] = obj.qp_solver.solve(obj.qpec);
            end

            % check is QP solve successful; use solution, or get a fallback (e.g. qpec solve)
            if stats.qp_solved
                % check is it better;
                if obj.mpc_options.use_feedback_qp
                    stats.qp_accepted = 1;
                    obj.feedback_phase_successful = true;
                    results_qpec = results_qp;
                elseif obj.mpc_options.use_probing_qp
                    f_k = full(obj.qpec_functions.f_fun(w_k, p_val));
                    f_k_new = full(obj.qpec_functions.f_fun(w_k+full(results_qp.x), p_val));
                    fprintf('QP objective %2.2f \n',f_k_new);
                    if f_k_new < obj.mpc_options.objective_ratio*f_k || norm(full(results_qp.x)) < 1e-6
                        stats.qp_accepted = 1;
                    else
                        stats.qp_accepted = 0;
                    end
                    results_qpec = results_qp;
                end
            else
                obj.feedback_phase_successful = false;
            end

            % If probing QP not good enough; solve QPCC
            if obj.mpc_options.use_probing_qp && ~stats.qp_accepted
                % In QP probing; if the QP is not good enough; solve the QPCC;
                fprintf('nosnoc mpc: Solving QPCC in feedback phase after QP probing was not accepted. \n')
                results_qpec = obj.qpec.solve();
                stats.qpec_solved = obj.qpec.stats.success;
                f_k_new = full(obj.qpec_functions.f_fun(w_k+full(results_qp.x), p_val));
                fprintf('QPCC objective %2.2f \n',f_k_new);
            end

            obj.ocp_solver.discrete_time_problem.f_result = f_k_new;

            % Check was the feedback phase computation a success
            if stats.qpec_solved || stats.qp_solved
                obj.feedback_phase_successful = true;
            else
                obj.feedback_phase_successful = false;
            end

            % New step; new multipliers
            if obj.feedback_phase_successful
                d_qpec    = full(results_qpec.x);
                lam_g_new = full(results_qpec.lam_g);
                lam_x_new = full(results_qpec.lam_x);
                % TODO: MPEC multipliers for G and H; xi and nu
                w_k       = w_k + d_qpec;

                if obj.mpc_options.warmstart_qpec
                    obj.qpec.w0 = d_qpec;
                    if isa(obj.qpec_options,'nosnoc.qpec.GurobiOptions')
                        obj.qpec.y0 = results_qpec.y; % y0 computed by updating qpec data! % only needed if we use gurobi
                    end
                end

                % note that multipliers are stored in ocp solver memory below
                stats.d_qpec_norm = norm(d_qpec);
                obj.stats.last_d_qpec_norm = stats.d_qpec_norm;
                % Update values in OCP solver memory
                obj.ocp_solver.discrete_time_problem.w.res = w_k;
                obj.ocp_solver.discrete_time_problem.w.mult = lam_x_new;
                obj.ocp_solver.discrete_time_problem.g.mult = lam_g_new;

                % Read first control and stats
                u_res = obj.ocp_solver.get('u');
                u = u_res(:,1);
                stats.ocp_solver_stats = obj.qpec.stats; % Why overwrite ocp solver stats with qpec stats? % Why store in ocp_solver_stats and in stats of SQPCC objec? -------> Mostly to initalize
            else
                % FALLBACK for FEEDBACK PHASE if QP/QPCC fails
                warning('nosnoc mpc: QPCC/QP step failed.');
                stats.ocp_solver_stats = struct; % return empty if failure happend;
                % keyboard;
                switch obj.mpc_options.on_qpec_failure
                    case 'resolve'
                        fprintf('nosnoc mpc: Resorting to full MPCC solve.\n');
                        obj.ocp_solver.set_x0(x0);
                        obj.ocp_solver.solve();
                        u_res = obj.ocp_solver.get('u');
                        u = u_res(:,1);
                        stats.d_qpec_norm = 0; % no QPCC step taken

                    case 'keep' % use_next: take the planned “next” control from last plan
                        fprintf('nosnoc mpc: Using next planned control (shifted) from previous OCP plan.\n');
                        u_plan = obj.ocp_solver.get('u');
                        if size(u_plan,2) >= 2
                            u = u_plan(:,2);         % the “next” control from previous plan
                        else
                            u = u_plan(:,1);         % fallback if horizon = 1
                        end
                        % Keep last known QPCC norm if you store it; otherwise set NaN
                        if isfield(obj.stats,'last_d_qpec_norm')
                            stats.d_qpec_norm = obj.stats.last_d_qpec_norm;
                        else
                            stats.d_qpec_norm = NaN;
                        end
                    case 'previous' % hold: re-apply last control actually used on the plant
                        fprintf('nosnoc mpc: Holding previous applied control.\n');
                        % You need to store last applied u after each get_feedback call:
                        % obj.stats.last_u_applied = u;  (set at the end of get_feedback)
                        if isfield(obj.stats,'last_u_applied')
                            u = obj.stats.last_u_applied;
                        else
                            % If not stored yet, fall back to first control of last plan
                            u_res = obj.ocp_solver.get('u');
                            u = u_res(:,1);
                        end
                        if isfield(obj.stats,'last_d_qpec_norm')
                            stats.d_qpec_norm = obj.stats.last_d_qpec_norm;
                        else
                            stats.d_qpec_norm = NaN;
                        end

                    otherwise
                        error('Invalid on_qpec_failure option: %s', obj.mpc_options.on_qpec_failure);
                end
                % end
            end
            stats.feedback_time = toc(feedback_timer);
            if isa(obj.qpec_options,'nosnoc.ccopt.Options')
                obj.stats.feedback_times{end+1} = obj.qpec.stats.ccopt.total_wall_time;
                stats.feedback_time = obj.qpec.stats.ccopt.total_wall_time;
                %stats.feedback_time = obj.qpec.stats.t_wall_total;
            else
                try
                    obj.stats.feedback_times{end+1} = obj.qpec.stats.wall_time_total;
                catch
                    % TODO; mpecopt may not return proper time. To be fixed.
                    obj.stats.feedback_times{end+1} = nan;
                end
                try
                    stats.feedback_time = obj.qpec.stats.wall_time_total;
                catch
                    % TODO; mpecopt may not return proper time. To be fixed.
                    stats.feedback_time = nan;
                end
            end

            if obj.feedback_phase_successful
                fprintf("nosnoc mpc: full_feedback_time: %2.5f, solver_time: %2.5f\n", stats.feedback_time, obj.stats.feedback_times{end})
            end

            obj.stats.last_u_applied = u; % if faild, can be reused

            if isfield(obj.stats, 'feedback_stats')
                obj.stats.feedback_stats(end+1) = stats;
            else
                obj.stats.feedback_stats = stats;
            end

        end

        function [x] = get_predicted_state(obj)
            % This method returns the predicted state from the last OCP solve. used for advanced step methods
            x_res = obj.ocp_solver.get("x");
            x = x_res(:,obj.problem_options.N_finite_elements(1)+1);
        end

        function [x] = get_current_state(obj)
            % This method returns the predicted state from the last OCP solve. used for advanced step methods
            x = obj.ocp_solver.discrete_time_problem.w.x(0,0,obj.ocp_solver.opts.n_s).lb;
            % x_res = obj.ocp_solver.get("x");
            % x = x_res(:,obj.problem_options.N_finite_elements(1));
        end

        % ----------------------------------- QPCC helper functions ---------------------------
        function qpec_functions = build_qpec_symbolic_functions(obj)
            % Build CasADi symbolic functions for QPCC evaluation
            % (gradients, Hessian, constraints) stored in the struct qpec_functions. (needed to set up a qpec);
            % Provides callable functions used to compute QPCC matrices numerically during MPC updates.
            import casadi.*
            obj.ocp_solver.solve(); % This does also the very first QPCC solver

            % get casadi variables and parameterse
            w = obj.ocp_solver.discrete_time_problem.w.sym;
            p = obj.ocp_solver.discrete_time_problem.p.sym;

            % get symbolic functions;
            % L_stage = o\bj.ocp_solver.discrete_time_problem.
            f = obj.ocp_solver.discrete_time_problem.f;
            g = obj.ocp_solver.discrete_time_problem.g.sym;
            G = obj.ocp_solver.discrete_time_problem.G.sym;
            H = obj.ocp_solver.discrete_time_problem.H.sym;

            % problem from mpecopt
            % w = obj.ocp_solver.discrete_time_problem.solver.mpec_casadi.x;
            % p = obj.ocp_solver.discrete_time_problem.solver.mpec_casadi.p;
            % f = obj.ocp_solver.discrete_time_problem.solver.mpec_casadi.f;
            % g = obj.ocp_solver.discrete_time_problem.solver.mpec_casadi.g;
            % G = obj.ocp_solver.discrete_time_problem.solver.mpec_casadi.x1;
            % H = obj.ocp_solver.discrete_time_problem.solver.mpec_casadi.x2;
            fun_opts.record_time = true;

            % Define functions for derivatives
            n_g = length(g);
            lam_g = SX.sym('lam_g', n_g);
            if obj.mpc_options.discard_constraints_in_hessian
                L = f;
            else
                L = f + lam_g'*g; % Lagrangian without linear terms (because we need only 2nd derivatives of L)
            end
            print_casadi_vector(L);

            % objective
            qpec_functions = struct;
            qpec_functions.f_fun = Function('f', {w, p}, {f}, fun_opts);
            qpec_functions.nabla_f_fun = Function('nabla_f', {w, p}, {gradient(f, w)}, fun_opts);
            qpec_functions.hess_L_fun = Function('hess_L', {w, lam_g, p}, {hessian(L, w)}, fun_opts);

            % regular constraints
            qpec_functions.g_fun = Function('g', {w, p}, {g}, fun_opts);
            qpec_functions.nabla_g_fun = Function('nabla_g', {w, p}, {jacobian(g, w)}, fun_opts);

            % complementarity constraints;
            qpec_functions.G_fun = Function('G', {w, p}, {G}, fun_opts);
            qpec_functions.H_fun = Function('H', {w, p}, {H}, fun_opts);
            qpec_functions.nabla_G_fun = Function('nabla_G', {w, p}, {jacobian(G, w)}, fun_opts);
            qpec_functions.nabla_H_fun = Function('nabla_H', {w, p}, {jacobian(H, w)}, fun_opts);
            qpec_functions.w = w;
        end

        function qpec = setup_qpec_solver(obj)
            % Setup the QPCC solver object using sparsity information and selected backend.
            % Allocates solver structure once and prepares it for numerical updates each iteration.
            % Backend is implicitly specified via obj.qpec_options
            import casadi.*
            % initalize QPCC (some solvers need sparsity patterns)
            qpec = nosnoc.qpec.Qpec(obj.qpec_functions.hess_L_fun.sparsity_out(0),...
                obj.qpec_functions.nabla_g_fun.sparsity_out(0),...
                obj.qpec_functions.nabla_G_fun.sparsity_out(0),...
                obj.qpec_functions.nabla_H_fun.sparsity_out(0));
            % create a qpec solver
            if isa(obj.qpec_options,'nosnoc.reg_homotopy.Options')
                qpec.create_qpec_solver(obj.qpec_options, 'reg_homotopy');
            elseif isa(obj.qpec_options,'nosnoc.qpec.GurobiOptions')
                % todo create gurobi solver
                qpec.create_qpec_solver(obj.qpec_options, 'gurobi');
            elseif isa(obj.qpec_options,'mpecopt.Options')
                qpec.create_qpec_solver(obj.qpec_options, 'mpecopt');
            elseif isa(obj.qpec_options,'nosnoc.ccopt.Options')
                qpec.create_qpec_solver(obj.qpec_options, 'ccopt');
            else
                % Otherwise it is a struct so assume it is lcqpow
                qpec.create_qpec_solver(obj.qpec_options, 'lcqpow');
            end
        end

        function setup_qp_solver(obj)
            % qp_solver_method = 'qp_ipopt', 'qp_highs' , could also be gurobi with appropate obj.qp_options (instead of directly via gurobi solver class)
            % Setup the QP solver using the Qpec class constructor
            % Allocates solver structure once and prepares it for numerical updates each iteration.
            if isa(obj.qp_options,'nosnoc.qpec.GurobiOptions')
                obj.qpec.create_qp_solver(obj.qp_options,'qp_gurobi');  % pass options and solverplug
            elseif isfield(obj.qp_options,'ipopt')
                obj.qpec.create_qp_solver(obj.qp_options,'qp_ipopt');  % pass options and solverplug
            elseif isfield(obj.qp_options,'highs')
                % check if highs
                % todo: implement highs interface
            elseif isempty(obj.qp_options)
                obj.qpec.create_qp_solver(obj.qp_options,'none'); % returns empty
            end

        end

        function initialize_qpec(obj)
            % do basic initalization of qpec;
            obj.qpec = obj.qpec;
            w_k = obj.ocp_solver.discrete_time_problem.w.res; % Read linearization point from ocp solver
            % w_k = obj.ocp_solver.discrete_time_problem.solver.solver_initialization.x0;
            % Bounds on general constriants
            % (no shift, as constraints read as lb <= g_k + nabla_g'*w <= ub)

            obj.qpec.lba = obj.ocp_solver.discrete_time_problem.g.lb;
            obj.qpec.uba = obj.ocp_solver.discrete_time_problem.g.ub;
            obj.qpec.w0 = zeros(size(w_k));  % initial guess for qpecs solver

            % default bounds on w
            obj.qpec.lbw = obj.ocp_solver.discrete_time_problem.w.lb;
            obj.qpec.ubw = obj.ocp_solver.discrete_time_problem.w.ub;

            % Eval G_k and H_k at sol and set up y0;
            p_val = obj.ocp_solver.discrete_time_problem.p.val;                     % problem parameters
            w_k = obj.ocp_solver.discrete_time_problem.w.res; % Read linearization point from ocp solver
            G_k = full(obj.qpec_functions.G_fun(w_k, p_val));
            H_k = full(obj.qpec_functions.H_fun(w_k, p_val));
            obj.qpec.y0 = G_k >= H_k;
        end

        function stats = update_qpec_data(obj)
            % Evaluate and update QPCC matrices (Q, A, G, H, q, b, g, h) from OCP.
            update_timer = tic;
            % obj.qpec = obj.qpec;
            p_val = obj.ocp_solver.discrete_time_problem.p.val;                     % problem parameters
            w_k = obj.ocp_solver.discrete_time_problem.w.res; % Read linearization point from ocp solver

            if obj.mpc_options.discard_constraints_in_hessian
                lam_g_k = zeros(size(obj.ocp_solver.discrete_time_problem.g.sym));
            else
                lam_g_k = obj.ocp_solver.discrete_time_problem.g.mult;
            end

            % Objective
            nabla_f_k = full(obj.qpec_functions.nabla_f_fun(w_k, p_val));
            Q_k = full(obj.qpec_functions.hess_L_fun(w_k, lam_g_k, p_val));

            % Regular constraints
            g_k = full(obj.qpec_functions.g_fun(w_k, p_val));
            nabla_g_k = full(obj.qpec_functions.nabla_g_fun(w_k, p_val));
            % Comp G
            G_k = full(obj.qpec_functions.G_fun(w_k, p_val));
            nabla_G_k = full(obj.qpec_functions.nabla_G_fun(w_k, p_val));
            % Comp H
            H_k = full(obj.qpec_functions.H_fun(w_k, p_val));
            nabla_H_k = full(obj.qpec_functions.nabla_H_fun(w_k, p_val));

            % update qpec data.
            Q_k = convexify_hessian(obj, Q_k);

            obj.qpec.Q = Q_k + obj.mpc_options.rho_d*eye(obj.qpec.dims.n_x);
            obj.qpec.A = nabla_g_k;
            obj.qpec.G = nabla_G_k;
            obj.qpec.H = nabla_H_k;
            obj.qpec.q = nabla_f_k;
            obj.qpec.b = g_k;
            obj.qpec.g = G_k;
            obj.qpec.h = H_k;
            % update y0; beacuse it was shifted!
            obj.qpec.y0 = G_k >= H_k;
            stats.update_time = toc(update_timer);
            % NOTE: This may be done in a generic QPCC, but here done in preparation phase after x0 is known
            % % Shift bounds  as lb <= w_k+d <= ub reads as lb - w_k <= d <= ub - w_k
            % obj.qpec_initialization.lbx = obj.ocp_solver.discrete_time_problem.w.lb - w_k;
            % obj.qpec_initialization.ubx = obj.ocp_solver.discrete_time_problem.w.ub - w_k;
        end

        function update_initial_value_in_qpec(obj)
            % Shift bounds  as lb <= w_k+d <= ub reads as lb - w_k <= d
            % <= ub - w_k; in MPC this is done after the new x0 is updated into w.lb, w.ub and w.res;
            obj.qpec.lbw = obj.ocp_solver.discrete_time_problem.w.lb - obj.ocp_solver.discrete_time_problem.w.res;
            obj.qpec.ubw = obj.ocp_solver.discrete_time_problem.w.ub - obj.ocp_solver.discrete_time_problem.w.res;
        end

        function Q_k = convexify_hessian(obj, Q_k)
            % Convexify Hessian matrix according to selected strategy.
            opts = obj.mpc_options;
            method = upper(opts.sqpec_hessian_convexification);

            switch method
                case 'NONE'
                    % nothing to do
                    return

                case {'PROJECT','MIRROR'}
                    % try Cholesky first
                    [~,p] = chol((Q_k+Q_k')/2);
                    if p == 0
                        % already PSD
                        return
                    end
                    % eigen-decomposition
                    [V,D] = eig(full((Q_k+Q_k')/2));
                    d = diag(D);
                    switch method
                        case 'PROJECT'
                            d = max(d, opts.eps_hessian);
                        case 'MIRROR'
                            d = abs(d);
                    end
                    Q_k = V*diag(d)*V';
                    % symmetrize for numerical safety
                    Q_k = 0.5*(Q_k+Q_k');

                case {'LEVENBERG_MARQUARDT'}
                    if isfield(opts,'rho_lm')
                        rho = opts.rho_lm;
                    else
                        rho = 1e1; % default mild regularization
                    end
                    n = size(Q_k,1);
                    Q_k = Q_k + rho*eye(n);

                case {'GERSHGORIN'}
                    Qsym = 0.5*(Q_k + Q_k'); % ensure symmetry
                    gersh_lb = min(diag(Qsym) - sum(abs(Qsym - diag(diag(Qsym))), 2));
                    if gersh_lb < opts.eps_hessian
                        rho = -gersh_lb + opts.eps_hessian;
                        Q_k = Qsym + rho*eye(size(Qsym));
                    else
                        Q_k = Qsym;
                    end

                otherwise
                    error('Unknown Hessian convexification method: %s', method);
            end
        end


        % --------------------------------- Getters and setters ----------------------------------

        function ret = get(obj,field)
            % Generic get wrapper for the ocp_solver underneath
            ret = obj.ocp_solver.get(field);
        end

        function ret = get_full(obj,field)
            % Generic get_full wrapper for the ocp_solver underneath
            ret = obj.ocp_solver.get_full(field);
        end

        function t_grid = get_time_grid(obj)
            t_grid = obj.ocp_solver.get_time_grid();
        end

        function t_grid_full = get_time_grid_full(obj)
            t_grid_full = obj.ocp_solver.get_time_grid_full();
        end

        function t_grid = get_control_grid(obj)
            t_grid = obj.ocp_solver.get_control_grid();
        end

        function set_param(obj, param, index, value)
            obj.ocp_solver.set_param(param, index, value);
        end

        function f_opt = get_objective(obj)
            % This method ocp objective
            f_opt = obj.ocp_solver.get_objective();
        end

    end
end
