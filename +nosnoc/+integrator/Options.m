classdef Options < handle

    properties
        %-------- FESD Integrator Specific -----------%
        
        fesd_solver_opts nosnoc.solver.Options % nosnoc FESD solver opts used for FESD solver. 
        use_previous_solution_as_initial_guess(1,1) logical = 0
        real_time_plot(1,1) logical = 0
        break_simulation_if_infeasible(1,1) logical = 0

        %--------- Smoothed PSS Integrator Specific -----------%
        
        matlab_ode_solver = 'ode23s';
        matlab_ode_opts = odeset;

        sigma_smoothing(1,1) double {mustBePositive} = 1e-7;
    end

    methods
        function obj = Options()

        end

        function [] = preprocess(obj)

        end

        function json = jsonencode(obj,varargin)
            import casadi.*
            opts_struct = struct(obj);

            json = jsonencode(opts_struct, varargin{:});
        end
    end
end
