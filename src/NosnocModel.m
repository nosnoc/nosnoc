classdef NosnocModel < handle

    properties
        %----- basic user input -----
        % state
        x
        x0
        lbx
        ubx

        % user algebraics
        z
        z0
        lbz
        ubz
        g_z % uzer algebraic constraints

        % global variables
        v_global
        v0_global
        lbv_global
        ubv_global

        % control
        u
        lbu
        ubu

        % global parameters
        p_global
        p_global_val

        % time varying parameters
        p_time_var
        p_time_var_val

        F         % Dynamic functions
        c         % Switching functions
        S         % Sign matrix
        g_Stewart % Stewart indicator functions
        
        f_q        % Stage cost
        f_q_T      % Terminal cost

        h % Step size
        h_k % Finite element step size

        N_finite_elements % Number of finite elements

        % least squares
        % TODO: is this necessary?

        % constraints
        g_path
        g_path_lb
        g_path_ub
        g_comp_path
        g_terminal

        %----- DCS/time_freezing mode user input -----
        f_c  % Gap functions
        mu_f % Friction coef
        e    % Restitution coef
        q    % Generalized position
        v    % Generalized velocity
        f_v  % Generalized forces
        M    % Inertia matrix

        J_normal
        J_tangent
        D_tangent

        %----- Algebraic variables -----

        % Stewart
        theta
        lam
        mu

        % Step
        alpha
        lambda_n
        lambda_p
        beta
        gamma

        % CLS
        lambda_normal
        y_gap
        lambda_tangent
        gamma_d
        beta_d
        delta_d
        beta_conic
        gamma_conic
        p_vt
        n_vt
        alpha_vt
        % CLS impulse vars
        Lambda_normal
        Y_gap
        Lambda_tangent
        Gamma_d
        Beta_d
        Delta_d
        Gamma_conic
        Beta_conic
        P_vt
        N_vt
        Alpha_vt

        %----- Generated model -----
        f_x % state dynamics (possibly not generated)

        g_switching
        g_z_all


        % Functions
        f_x_fun
        g_z_all_fun
        g_switching_fun
        c_fun
        dot_c_fun
        g_Stewart_fun
        g_impulse_fun
        f_c_fun
        M_fun
        invM_fun
        J_normal_fun
        J_tangent_fun
        D_tangent_fun
        f_q_T_fun
        J_cc_fun
        % TODO lsq

        % params
        p_time_var_stages
        p_dyn

        % Dimensions
        dims
    end

    methods
        function obj = NosnocModel()

        end

        function vars = generate_vars(obj,settings)

        end

        function reformulate(obj, settings)

        end

        function verify_and_backfill(obj, settings)
            if isfield(obj,'x')

            end
        end
        
    end

end
