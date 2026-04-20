% This file is part of nosnoc.
classdef Options < handle
    % Options for the Sequential QPEC (SQPEC) algorithm.
    % Stores configuration parameters controlling numerical tolerances,
    % globalization, verbosity, and solver-specific options.

    properties
        % ---- General ----
        solver_name {mustBeTextScalar} = 'sqpec_solver'                 % Name identifier for this solver instance.
        casadi_symbolic_mode {mustBeMember(casadi_symbolic_mode,{'casadi.SX', 'casadi.MX'})} = 'casadi.SX'  % Symbolic type used for CasADi function generation.

        % ---- MPCC / QPEC and Homotopy settings ----
        complementarity_tol (1,1) double {mustBeReal, mustBePositive} = 1e-9    % Tolerance for complementarity satisfaction.
        tol (1,1) double {mustBeReal, mustBePositive} = 1e-8                    % General numerical tolerance for convergence checks.
        assume_lower_bounds (1,1) logical = true                                % Assume all inequality constraints have lower bounds (optimization shortcut).
        lift_complementarities (1,1) logical = false                            % Whether to lift complementarity constraints into auxiliary variables.
        initial_comp_all_zero (1,1) logical = false                             % Initialize complementarity variables to zero (for debugging or special tests).
        discard_constraints_in_hessian (1,1) logical = true                     % Ignore constraint contributions in Hessian computation.
        max_iter (1,1) double {mustBeInteger, mustBePositive} = 20              % Maximum number of SQPEC iterations.
        qpec_solver_options = nosnoc.reg_homotopy.Options()                     % Options object passed to the internal QPEC solver (e.g., reg_homotopy, mpecopt).
        rho_d (1,1) double {mustBeReal, mustBeNonnegative} = 0                  % Regularization weight added to the QPEC Hessian (ρ_d * I).

        % ---- Verbosity ----
        print_level = 3                                                         % Print verbosity level (0 = silent, higher = more detail).
        print_details_if_infeasible = 0                                         % Print additional diagnostics if infeasibility is detected.
        pause_homotopy_solver_if_infeasible = 0                                 % Pause execution if the homotopy solver fails (debugging).

        % ---- Globalization / Filter line search ----
        use_globalization (1,1) logical = true                                  % Enable globalization via filter line search.
        filter_gamma_f (1,1) {mustBeInRange(filter_gamma_f,0,1,'exclusive')} = 1e-5   % Filter parameter γ_f for objective reduction.
        filter_gamma_h (1,1) {mustBeInRange(filter_gamma_h,0,1,'exclusive')} = 1e-5   % Filter parameter γ_h for constraint violation.
        filter_s_f (1,1) {mustBeInRange(filter_s_f,1,50)} = 2.3                      % Filter scaling parameter s_f (objective).
        filter_s_h (1,1) {mustBeInRange(filter_s_h,1,51,'exclusive')} = 1.1          % Filter scaling parameter s_h (constraints).
        filter_eta_f (1,1) {mustBeInRange(filter_eta_f,0,0.5,'exclusive')} = 1e-4    % Filter acceptance parameter η_f ∈ (0,0.5).
        filter_h_min (1,1) double {mustBeReal, mustBeNonnegative} = 1e-7             % Minimum constraint violation threshold for filter.
        filter_add_all_acceptable_points (1,1) logical = false                       % If true, add all acceptable points to the filter (otherwise only strict improvements).

        % ---- Timeouts (optional, currently unused) ----
        % timeout_cpu (1,1) {mustBeReal, mustBeNonnegative} = 0                      % CPU time limit [s].
        % timeout_wall (1,1) {mustBeReal, mustBeNonnegative} = 0                     % Wall-clock time limit [s].
    end
end
