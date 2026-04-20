classdef BnlpLpcc < nosnoc.mpc.Base
    % Armin: This does not work well - we should probably delete it!
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
                % obj.solver_options.initialization_strategy = 'TakeInitialGuessActiveSet';
                if 1
                    obj.solver_options.initialization_strategy = 'TakeInitialGuessDirectly';
                    obj.solver_options.initialization_strategy = 'TakeInitialGuessActiveSet';

                    w_opt = obj.ocp_solver.get_w();
                    w_opt = obj.ocp_solver.get_control_grid();
                    % theta_opt = obj.ocp_solver.get('theta');
                    % lambda_opt = obj.ocp_solver.get('lambda');
                    x = obj.ocp_solver.get('x');
                    % rho_guess = 2*(max(max(theta_opt(:),lambda_opt(:)))+1);
                    % % rho_guess = max(abs(w_opt));
                    rho_guess = 2*norm(x(:,1)-x0,inf)+1;
                    % rho_guess = norm(x(:,1)-x0,inf)+1;
                    % obj.solver_options.rho_TR_phase_ii_init = rho_guess;
                    obj.solver_options.rho_TR_phase_ii_init = min(0.5,rho_guess);
                    obj.solver_options.max_iter = 1;
                    obj.solver_options.max_inner_iter = 1;
                else
                    obj.solver_options.relax_and_project_homotopy_parameter_steering = "Ell_inf";
                    obj.solver_options.rho_TR_phase_ii_init = 1e-3;
                    obj.solver_options.rho_TR_phase_i_init = 1e-3;
                    obj.solver_options.relax_and_project_sigma0 = 1E-3;
                    obj.solver_options.initialization_strategy = 'RelaxAndProject';
                end
            else
                obj.solver_options.initialization_strategy = 'RelaxAndProject';
                % obj.solver_options.rho_TR_phase_ii_init = 1e-2;
                % obj.solver_options.max_iter = 5;
                % obj.solver_options.max_inner_iter = 6;
            end
            obj.ocp_solver.set_x0(x0);
            obj.ocp_solver.solve();
            u_res = obj.ocp_solver.get('u');
            u = u_res(:,1);
            stats.ocp_solver_stats = obj.ocp_solver.stats;
            stats.feedback_time = toc(feedback_timer);
            obj.last_solve_successful = true;
        end

        function [stats] = do_preparation(obj)
            % This method should be called after get feedback to do the work of
            % setting up the warm starting, it does nothing if the previous solve
            % did not converge.
            preparation_timer = tic;
            if obj.ocp_solver.stats.h_total < 1e-8 % TODO: param here
                if obj.mpc_options.do_shift_initialization
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
