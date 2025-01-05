classdef Options < handle

    properties
        integrator_plugin IntegratorType = IntegratorType.FESD; % Type of integrator to use.
        
        %-------- FESD Integrator Specific -----------%
        
        fesd_solver_opts nosnoc.solver.Options = nosnoc.solver.Options() % nosnoc FESD solver opts used for FESD solver. 
        use_previous_solution_as_initial_guess(1,1) logical = 0
        real_time_plot(1,1) logical = 0
        break_simulation_if_infeasible(1,1) logical = 0

        %--------- Smoothed PSS Integrator Specific -----------%
        
        matlab_ode_solver = 'ode23s';
        matlab_ode_opts = odeset;

        sigma_smoothing(1,1) double {mustBePositive} = 1e-7;

        %--------- Stewart's PSS Integrator Specific -----------%
        alpha_stewart (1,1) double {mustBePositive} = 1e1; % positive constant for ensuring that the matrix M is copositive.
        epsilon_stewart (1,1) double {mustBePositive} = 1e-8; % epsilon for determining the epsilon active set in Stewart's method.
    end

    methods
        function obj = Options()
        end

        function [] = preprocess(obj)
            if obj.integrator_plugin == IntegratorType.SMOOTHED_PSS && ~ismember(obj.matlab_ode_solver, {'ode23s', 'ode15s', 'cvodesstiff', 'idas'})
                nosnoc.warning('likely_bad_ode_solver', 'You have selected an ode solver that is unlikely to be efficient for PSS due to stiffness. Choosing ode23s or ode15s are likely to work better. In MATLAB >=2024a cvodesstiff or idas can be more efficient than ode23s and ode15s.')
            end
        end

        function json = jsonencode(obj,varargin)
            import casadi.*
            opts_struct = struct(obj);

            json = jsonencode(opts_struct, varargin{:});
        end
    end
end
