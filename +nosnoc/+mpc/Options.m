classdef Options < handle
    properties
        fast_sigma_0(1,1) double {mustBeNonnegative} = 1e0;
        do_shift_initialization(1,1) logical = true;
        rho_d(1,1) double {mustBeNonnegative} = 0; % Penalty parameter in rho_d*d'*d, to keep solution closer to previous. essentialy a hessian regularization
        discard_constraints_in_hessian (1,1) logical = true; % If True, the Lagrangian Hessian is equal to nabla^2 f(w). If linear least square objective, then equal to Gauss-Newton.

        % advanced mpc settings
        solve_advanced_problem (1,1) logical = false; % In preparation phase solve a problem with a predicted state;
        advanced_problem_type {mustBeMember(advanced_problem_type,{'full','sqpec'})} = 'sqpec';
        advanced_n_qpecs  (1,1) double {mustBeInteger, mustBePositive} = 3 

        on_qpec_failure = 'previous'; % 'resolve' | 'keep' | 'previous'% keep   -> use next planned control (u(:,2)) % previous -> hold last applied control

        sqpec_hessian_convexification (1,1) nosnoc.sqpec.HessianConvexification = 'NONE' % What method to convexify the QPEC hessian.
        rho_lm(1,1) double {mustBeNonnegative} = 1e1;
        eps_hessian(1,1) double {mustBeNonnegative} = 1e-6;
        warmstart_qpec(1,1) logical = true;
        warmstart_full_mpc(1,1) logical = true;
        warmstart_mid_homotopy(1,1) logical = true; % Used?!
        objective_ratio(1,1) double {mustBeNonnegative} = 0.97; % value of new obj as fraction of old for qp probing to be accepted;

        use_probing_qp(1,1) logical = false; % Attemp to solve qp with latest avilable active set; if fails fall back to qpcc;
        use_feedback_qp(1,1) logical = false; % In advanced step algortihms; solve only QP in the feedback;
    end
end
