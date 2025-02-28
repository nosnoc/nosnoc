classdef BnlpLpcc < nosnoc.mpc.Base
% An MPC implementation using a single outer iteration of the mpecopt solver.
% I.e. at each iteration it solves a single lpec and BNLP
    properties (Access=private)
        ocp_solver nosnoc.ocp.Solver % nosnoc OCP solver used to generate feedback.
        mpc_options nosnoc.mpc.Options % nosnoc mpc specific options.
        problem_options nosnoc.Options % nosnoc problem options.
        solver_options mpecopt.Options % mpecopt solver options.
        model % nosnoc model.
        
        last_solve_successful(1,1) logical = false % true if last solve was successful, false otherwise. 
    end
    
    methods
        function obj = BnlpLpcc(model, mpc_options, problem_options, solver_options)
            obj.mpc_options = mpc_options;
            obj.problem_options = problem_options;
            obj.solver_options = solver_options;
            obj.ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
        end
        
        function [u, stats] = get_feedback(obj, x0)
            % This method takes a state estimate $x_0$ and returns the corresponding control $u$.
            feedback_timer = tic;

            % Update sigma_0 to the fast one if we have 
            if obj.last_solve_successful
                solver_options.initialization_strategy = 'TakeInitialGuessDirectly';
                solver_options.max_iter = 1;
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
            if isfield(obj.ocp_solver.stats, "converged") && obj.ocp_solver.stats.converged
                obj.ocp_solver.do_warmstart();
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
