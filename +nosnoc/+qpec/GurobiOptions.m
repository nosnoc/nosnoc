% This file is part of nosnoc.
classdef GurobiOptions < handle
    % Options for Gurobi QPCC solver (MIQP, SOS1, or Scholtes regularization)
    %
    % These parameters control the reformulation type and solver behavior.
    % Native Gurobi parameters can be passed via `gurobi_params`.

    properties
        % --- High-level method selection ---
        method {mustBeMember(method,{'miqp','bigm','sos1','reg','qp'})} = 'miqp';
        % Type of QPCC reformulation used:
        % 'miqp'  -> Big-M mixed-integer QP
        % 'sos1'  -> Complementarity as SOS1 sets
        % 'reg'   -> Scholtes regularization (continuous relaxation)
        tau0 (1,1) double {mustBeReal, mustBePositive} = 1.0                                  % Initial relaxation τ.
        kappa (1,1) double {mustBeReal, mustBePositive} = 0.1                               % Homotopy reduction factor.
        N_homotopy (1,1) double {mustBeInteger, mustBePositive} = 9                          % # outer homotopy stages.
        Nsqp (1,1) double {mustBeInteger, mustBePositive} = 10                               % # SQP iterations per stage.
        sqp_tol (1,1) double {mustBeReal, mustBePositive} = 1e-8                             % SQP convergence tolerance.
        verbose (1,1) logical = true                                                         % Print per-stage progress.
        warmstart (1,1) logical = true                                                       % Enable warm-starting if available.
        warmstart_qp(1,1) logical = false                                                    % when solving qp with fixed y, warm start with avilable qpec.w0
        % miqp
        bigM (1,1) double {mustBePositive} = 1e2;                                            % Big-M constant for MIQP
        recover_duals (1,1) logical = true                                                   % Read duals from gurobi sol
        resolve_qp_for_duals (1,1) logical = true                                           % If true, resolve qp for fixed y to get duals, essentially resolve branch QP. If false, take MIQP duals.
        lower_bounds_comps (1,1) logical = false                                            % Add lower bounds on comps; or asume they aredone elserwhere (in simple bounds, a fesd thing)


        % --- Native Gurobi parameters ---
        gurobi_params = struct( ...
            'OutputFlag', 0, ...       % no console output
            'MIPGap', 1e-8, ...        % relative MIP optimality gap
            'Presolve', 0, ...         % full presolve % 2 is argessive, - 1 is auto
            'Method', 3, ...           % use Barrier method for continuous relaxations % 2 is barrier ; 4 is coccurent
            'NodeLimit', 2e3, ...       % limit number of explored nodes
            'FeasibilityTol',1e-2,...
            'NumericFocus',2 ...
            );
            % 'IntFeasTol',1e-2...

        % These are passed directly to gurobi(model, gurobi_params)
    end
end
