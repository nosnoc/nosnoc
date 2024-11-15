classdef FullMpcc < nosnoc.mpc.Base
% A basic MPC implementation which uses a call to the full MPCC solver to generate feedback.
% It does either shift or in place warmstarting during the preparation phase and uses a possibly reduced
% initial relaxation in the homotopy solver if the prior solution converged.
    properties (Access=private)
        ocp_solver nosnoc.ocp.Solver % nosnoc OCP solver used to generate feedback.
        mpc_options nosnoc.mpc.Options % nosnoc mpc specific options.
        problem_options nosnoc.Options % nosnoc problem options.
        solver_options nosnoc.solver.Options % nosnoc solver options.
        model % nosnoc model.

        last_solve_successful(1,1) logical = false % true if last solve was successful, false otherwise.
        cold_sigma_0(1,1) double {mustBeNonnegative} % The cold start sigma from solver_options. 
    end
    
    methods
        function obj = FullMpcc(model, mpc_options, problem_options, solver_options)
            obj.mpc_options = mpc_options;
            obj.problem_options = problem_options;
            obj.solver_options = solver_options;
            obj.cold_sigma_0 = solver_options.sigma_0;
            obj.ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);

            if mpc_options.fullmpcc_progressive_relaxation ~= 0
                obj.solver_options.progressive_relaxation_factor = mpc_options.fullmpcc_progressive_relaxation;
            end
        end
        
        function [u, stats] = get_feedback(obj, x0)
            % This method takes a state estimate $x_0$ and returns the corresponding control $u$.
            feedback_timer = tic;

            % Update sigma_0 to the fast one if we have 
            if obj.last_solve_successful
                obj.solver_options.sigma_0 = obj.mpc_options.fullmpcc_fast_sigma_0;
                obj.solver_options.N_homotopy = obj.mpc_options.fullmpcc_fast_N_homotopy;
            end
            obj.ocp_solver.set_x0(x0);
            obj.ocp_solver.solve();
            u_res = obj.ocp_solver.get('u');
            u = u_res(:,1);
            stats.ocp_solver_stats = obj.ocp_solver.stats;
            stats.feedback_time = toc(feedback_timer);
        end

        function [stats] = do_preparation(obj)
            % This method should be called after get feedback to do the work of
            % setting up the warm starting, it does nothing if the previous solve
            % did not converge.
            preparation_timer = tic;
            if true%isfield(obj.ocp_solver.stats, "converged") && obj.ocp_solver.stats.converged
                if obj.mpc_options.fullmpcc_do_shift_initialization
                    obj.ocp_solver.do_shift_initialization();
                else
                    obj.ocp_solver.do_warmstart();
                end
                obj.last_solve_successful = true;
            else
                obj.last_solve_successful = false;
            end
            stats.preparation_time = toc(preparation_timer);
        end

        function [x] = get_predicted_state(obj)
            % This method returns the predicted state from the last solve.
            x_res = obj.ocp_solver.get("x");
            x = x_res(:,obj.problem_options.N_finite_elements(1)+1);
        end

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
    end
end

