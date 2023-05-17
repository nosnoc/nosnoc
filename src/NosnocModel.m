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
        u0
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
        lsq_x
        lsq_u
        lsq_T

        % constraints
        g_path
        g_path_lb
        g_path_ub
        g_comp_path
        g_terminal

        % Terminal time
        T

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

        % all
        z_all
        
        % Stewart
        theta
        theta_sys
        lambda
        lambda_sys
        mu
        mu_sys

        % Step
        alpha
        alpha_sys
        lambda_n
        lambda_n_sys
        lambda_p
        lambda_p_sys
        beta
        gamma
        theta_step      % casadi expression either lifted or not
        theta_step_sys

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
        f_lsq_x_fun
        x_ref_val
        f_lsq_u_fun
        u_ref_val
        f_lsq_T_fun
        x_ref_end_val
        
        % params
        p_time_var_stages
        p_dyn

        % flags
        friction_exists
        
        % Dimensions
        dims
    end

    methods
        function obj = NosnocModel()
            obj.dims = NosnocDimensions();
        end

        function reformulate(obj, settings)

        end

        function vars = generate_vars(obj,settings)
            import casadi.*
            casadi_symbolic_mode = settings.casadi_symbolic_mode;
            dcs_mode = settings.dcs_mode;
            dims = obj.dims;
            g_lift_theta_step = [];
            g_lift_beta = [];
            switch dcs_mode
              case 'Stewart'
                % dimensions
                n_theta = sum(obj.dims.n_c_sys); % number of modes
                n_lambda = n_theta;
                for ii = 1:dims.n_sys
                    ii_str = num2str(ii);
                    % define theta (Filippov multiplers)
                    obj.theta_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['theta_' ii_str],obj.dims.n_c_sys(ii));
                    obj.theta = [obj.theta;obj.theta_sys{ii}];
                    % define mu_i (Lagrange multipler of e'theta =1;)
                    obj.mu_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['mu_' ii_str],1);
                    obj.mu = [obj.mu;obj.mu_sys{ii}];
                    % define lambda_i (Lagrange multipler of theta >= 0;)
                    obj.lambda_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['lambda_' ii_str],obj.dims.n_c_sys(ii));
                    obj.lambda = [obj.lambda;obj.lambda_sys{ii}];
                end
              case 'Step'
                n_alpha = sum(obj.dims.n_c_sys);
                n_lambda_n = sum(obj.dims.n_c_sys);
                n_lambda_p = sum(obj.dims.n_c_sys);
                % for creae_nlp_fesd
                n_theta = 2*n_alpha;
                n_lambda = n_lambda_n+n_lambda_p;
        
                for ii = 1:dims.n_sys
                    ii_str = num2str(ii);
                    % define alpha (selection of a set valued step function)
                    if ~settings.general_inclusion
                        obj.alpha_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['alpha_' ii_str],obj.dims.n_c_sys(ii));
                        obj.alpha = [obj.alpha;obj.alpha_sys{ii}];
                    else
                        % TODO this needs to change if subsystems.
                        obj.alpha_sys{ii} = obj.alpha{ii};
                    end
                    % define lambda_n_i (Lagrange multipler of alpha >= 0;)
                    obj.lambda_n_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['lambda_n_' ii_str],obj.dims.n_c_sys(ii));
                    obj.lambda_n = [obj.lambda_n;obj.lambda_n_sys{ii}];
                    % define lambda_p_i (Lagrange multipler of alpha <= 1;)
                    obj.lambda_p_sys{ii} = define_casadi_symbolic(casadi_symbolic_mode,['lambda_p_' ii_str],obj.dims.n_c_sys(ii));
                    obj.lambda_p = [obj.lambda_p;obj.lambda_p_sys{ii}];
                end

                % Define already here lifting variables and functions
                % TODO allow for custom beta lifting
                % Theta collects the vector for dot_x = F(x)Theta,
                % terms or theta_step from lifting;
                if ~settings.general_inclusion
                    for ii = 1:dims.n_sys
                        theta_temp = [];
                        ii_str = num2str(ii);
                        S_temp = obj.S{ii};
                        if settings.pss_lift_step_functions
                            % TODO implement automatic lifting
                        else
                            if ~settings.time_freezing_inelastic
                                for j = 1:size(S_temp,1)
                                    alpha_ij = 1;
                                    for k = 1:size(S_temp,2)
                                        % create multiaffine term
                                        if S_temp(j,k) ~=0
                                            alpha_ij = alpha_ij * (0.5*(1-S_temp(j,k))+S_temp(j,k)*obj.alpha_sys{ii}(k) ) ;
                                        end
                                    end
                                    theta_temp = [theta_temp;alpha_ij];
                                end
                            end
                        end
                        obj.theta_step_sys{ii} = theta_temp;
                    end
                end

                %% time-freezing inelastic impacts (exploit structure with taiolored formulae)
                if settings.time_freezing_inelastic
                    % theta_step are the lifting variables that enter the ODE r.h.s.
                    if ~settings.nonsmooth_switching_fun
                        alpha_q = alpha(1:dims.n_contacts);
                        alpha_v_normal = alpha(dims.n_contacts+1:2*dims.n_contacts);
                        if friction_exists
                            alpha_v_tangent = obj.alpha(2*dims.n_contacts+1:end);
                        end
                    else
                        alpha_qv = alpha(1:dims.n_contacts);
                        if friction_exists
                            alpha_v_tangent = obj.alpha(dims.n_contacts+1:end);
                        end
                    end

                    obj.theta_step = define_casadi_symbolic(casadi_symbolic_mode,'theta_step',n_aux+1);
                    theta_step_expr = eval([casadi_symbolic_mode '.zeros(' num2str(n_aux) '+1,1);']); % expressions for Filippov multipliers via alpha (and possibley beta).

                    % lifting variables
                    % empty expressions for initalization
                    beta_bilinear_ode = []; % for lifting bilinear terms in free flight dynamics multiplire
                    beta_bilinear_aux = []; % for lifiting bilinear terms appearing in aux. dynamics mutiplieres
                    beta_prod = []; % for lifting the multi affine term defineng the overall free flight dynamics multpliers
                                    % expressions for lifting
                    beta_bilinear_ode_expr = [];
                    beta_bilinear_aux_expr = [];
                    beta_prod_expr = [];
                    beta_prod_expr_guess = []; % extra expresion to make depend only on alpha (the one above depens on both and alpha and beta) - needed for eval. of inital guess

                    if pss_lift_step_functions
                        % lift bilinear terms in product terms for free flight ode % (alpha_q*alpha_v)
                        if ~nonsmooth_switching_fun
                            beta_bilinear_ode = define_casadi_symbolic(casadi_symbolic_mode,'beta_bilinear_ode',dims.n_contacts);
                            beta_bilinear_ode_expr = eval([casadi_symbolic_mode '.zeros(' num2str(dims.n_contacts) ',1);']);
                            if friction_exists
                                % lift bilinear terms defining aux dynamics (1-alpha_q)*(1-alpha_v)
                                beta_bilinear_aux = define_casadi_symbolic(casadi_symbolic_mode,'beta_bilinear_aux',dims.n_contacts);
                                beta_bilinear_aux_expr = eval([casadi_symbolic_mode '.zeros(' num2str(dims.n_contacts) ',1);']);
                            end
                        end
                        if dims.n_contacts > 2
                            beta_prod = define_casadi_symbolic(casadi_symbolic_mode,'beta_bilinear',dims.n_contacts-2);
                            beta_prod_expr = eval([casadi_symbolic_mode '.zeros(' num2str(dims.n_contacts) '-2,1);']);
                            beta_prod_expr_guess = eval([casadi_symbolic_mode '.zeros(' num2str(dims.n_contacts) '-2,1);']);
                        end
                    end
                    obj.beta = [beta_bilinear_ode;
                        beta_bilinear_aux;
                        beta_prod];
                    % expresions for theta's and lifting
                    %% Filippov multipliers
                    alpha_ode = 1; % initalized product for free flight multiplier
                    if ~pss_lift_step_functions
                        for ii = 1:dims.n_contacts
                            if nonsmooth_switching_fun
                                alpha_ode = alpha_ode*alpha_qv(ii);
                                if friction_exists
                                    theta_step_expr(ii+1) = (1-alpha_qv(ii))*(alpha_v_tangent(ii));
                                    theta_step_expr(1+dims.n_contacts+ii) = (1-alpha_qv(ii))*(1-alpha_v_tangent(ii));
                                else
                                    theta_step_expr(ii+1)=(1-alpha_qv(ii));
                                end
                            else
                                alpha_ode = alpha_ode*(alpha_q(ii)+alpha_v_normal(ii)-alpha_q(ii)*alpha_v_normal(ii));
                                if friction_exists
                                    theta_step_expr(ii+1) = (1-alpha_q(ii))*(1-alpha_v_normal(ii))*(alpha_v_tangent(ii));
                                    theta_step_expr(1+dims.n_contacts+ii) = (1-alpha_q(ii))*(1-alpha_v_normal(ii))*(1-alpha_v_tangent(ii));
                                else
                                    theta_step_expr(ii+1) = (1-alpha_q(ii))*(1-alpha_v_normal(ii));
                                end
                            end
                        end
                        theta_step_expr(1) = alpha_ode;
                    else
                        % lift and have bilinear terms
                        if nonsmooth_switching_fun
                            if dims.n_contacts <= 2
                                for ii = 1:dims.n_contacts
                                    alpha_ode = alpha_ode*alpha_qv(ii);
                                end
                                theta_step_expr(1) = alpha_ode;
                            else
                                beta_prod_expr(1) = (alpha_qv(1))*(alpha_qv(2));
                                beta_prod_expr_guess(1) = (alpha_qv(1))*(alpha_qv(2));
                                % lifting terms in between
                                for ii = 3:dims.n_contacts-1
                                    beta_prod_expr(ii-1) = beta_prod(ii-2)*(alpha_qv(ii)); % beta_{i} = beta{i-1}*(prod_term_i+1}
                                    beta_prod_expr_guess(ii-1) = beta_prod_expr_guess(ii-2)*(alpha_qv(ii)); % this is to have an expression depending only on alpha for the inital guess eval
                                end
                                theta_step_expr(1)= beta_prod(end)*(alpha_qv(dims.n_contacts)); % last lifting term;
                            end
                            % lifting of aux dyn multiplier expressions
                            for ii = 1:dims.n_contacts
                                if friction_exists
                                    theta_step_expr(ii+1) = (1-alpha_qv(ii))*(alpha_v_tangent(ii));
                                    theta_step_expr(1+dims.n_contacts+ii) = (1-alpha_qv(ii))*(1-alpha_v_tangent(ii));
                                else
                                    theta_step_expr(ii+1)=(1-alpha_qv(ii));
                                end
                            end
                        else
                            % with two smooth switching functions
                            beta_bilinear_ode_expr = alpha_q.*alpha_v_normal;
                            if friction_exists
                                beta_bilinear_aux_expr = (1-alpha_q).*(1-alpha_v_normal);
                            end

                            if dims.n_contacts <= 2
                                % here no lifting of product terms
                                for ii = 1:dims.n_contacts
                                    alpha_ode = alpha_ode*(alpha_q(ii)+alpha_v_normal(ii)-beta_bilinear_ode(ii));
                                end
                                theta_step_expr(1) = alpha_ode;
                            else
                                % here lifting of product terms
                                g_z_tf_beta_prod  = [beta_prod(1) - (alpha_q(1)+alpha_v_normal(1)-beta_bilinear_ode(1))*(alpha_q(2)+alpha_v_normal(2)-beta_bilinear_ode(2))]; % first lifting terms
                                                                                                                                                                              % lifting terms in between
                                for ii = 3:dims.n_contacts-1
                                    % beta_{i} = beta{i-1}*(prod_term_i+1}
                                    beta_prod_expr(ii-1) = beta_prod(ii-2)*(alpha_q(ii)+alpha_v_normal(ii)-beta_bilinear_ode(ii));
                                    beta_prod_expr_guess(ii-1) = beta_prod_expr_guess(ii-2)*(alpha_q(ii)+alpha_v_normal(ii)-beta_bilinear_ode(ii));
                                end
                                % last lifting term;
                                theta_step_expr(1) = beta_prod(end)*(alpha_q(dims.n_contacts)+alpha_v_normal(dims.n_contacts)-beta_bilinear_ode(dims.n_contacts));
                            end
                            % lifting of aux dyn multiplier expressions
                            for ii = 1:dims.n_contacts
                                if friction_exists
                                    theta_step_expr(ii+1) =  beta_bilinear_aux(ii)*(alpha_v_tangent(ii));
                                    theta_step_expr(1+dims.n_contacts+ii) = beta_bilinear_aux(ii)*(1-alpha_v_tangent(ii));
                                else
                                    theta_step_expr(ii+1) = (1-alpha_q(ii))*(1-alpha_v_normal(ii));
                                end
                            end
                        end
                    end
                    % equality constraints in DCS
                    g_lift_theta_step = obj.theta_step-theta_step_expr;
                    g_lift_beta = beta - [beta_bilinear_ode_expr; beta_bilinear_aux_expr; beta_prod_expr];
                    % auxiliary functions to get inital guess for new algebraic variables theta and beta
                    g_lift_theta_step_fun  = Function('g_lift_theta_step_fun',{alpha,beta},{theta_step_expr});
                    g_lift_beta_fun = Function('g_lift_beta_fun',{alpha},{[beta_bilinear_ode_expr;beta_bilinear_aux_expr;beta_prod_expr_guess]});
                    obj.theta_step_sys{1} = theta_step;
                end
                n_beta = length(obj.beta);
                n_theta_step = length(obj.theta_step);
              case 'CLS'
                lambda_normal = define_casadi_symbolic(casadi_symbolic_mode,'lambda_normal',n_contacts);
                y_gap = define_casadi_symbolic(casadi_symbolic_mode,'y_gap',n_contacts);
                % Variables for impulse equations
                Lambda_normal = define_casadi_symbolic(casadi_symbolic_mode,'Lambda_normal',n_contacts);
                Y_gap = define_casadi_symbolic(casadi_symbolic_mode,'Y_gap',n_contacts);
                %         P_vn = define_casadi_symbolic(casadi_symbolic_mode,'P_vn',n_contacts); % pos part of state jump law
                %         N_vn = define_casadi_symbolic(casadi_symbolic_mode,'N_vn',n_contacts); % neg part of state jump law
                L_vn = define_casadi_symbolic(casadi_symbolic_mode,'L_vn',n_contacts); % lifting variable for state jump law
                if friction_exists
                    % tangetial contact froce (firction force)
                    lambda_tangent = define_casadi_symbolic(casadi_symbolic_mode,'lambda_tangent',n_tangents);
                    % Impulse varaibles
                    Lambda_tangent = define_casadi_symbolic(casadi_symbolic_mode,'Lambda_tangent',n_tangents);
                    if isequal(friction_model,'Polyhedral')
                        gamma_d = define_casadi_symbolic(casadi_symbolic_mode,'gamma_d',n_contacts);
                        beta_d = define_casadi_symbolic(casadi_symbolic_mode,'beta_d',n_contacts); % lift friction cone bound
                        delta_d = define_casadi_symbolic(casadi_symbolic_mode,'delta_d',n_tangents); % lift lagrangian
                                                                                                     % Impulse varaibles
                        Gamma_d = define_casadi_symbolic(casadi_symbolic_mode,'Gamma_d',n_contacts);
                        Beta_d = define_casadi_symbolic(casadi_symbolic_mode,'Beta_d',n_contacts); % lift friction cone bound
                        Delta_d = define_casadi_symbolic(casadi_symbolic_mode,'Delta_d',n_tangents); % lift lagrangian
                    end
                    if isequal(friction_model,'Conic')
                        gamma = define_casadi_symbolic(casadi_symbolic_mode,'gamma',n_contacts);
                        beta = define_casadi_symbolic(casadi_symbolic_mode,'beta',n_contacts);
                        % Impulse variables;
                        Gamma = define_casadi_symbolic(casadi_symbolic_mode,'Gamma',n_contacts);
                        Beta = define_casadi_symbolic(casadi_symbolic_mode,'Beta',n_contacts);
                        switch conic_model_switch_handling
                          case 'Plain'
                            % no extra constraints
                          case 'Abs'
                            p_vt = define_casadi_symbolic(casadi_symbolic_mode,'p_vt',n_tangents); % positive parts of tagnetial velocity (for switch detection)
                            n_vt = define_casadi_symbolic(casadi_symbolic_mode,'n_vt',n_tangents); % negative parts of tagnetial velocity (for switch detection)
                                                                                                   % Impulse
                            P_vt = define_casadi_symbolic(casadi_symbolic_mode,'P_vt',n_tangents);
                            N_vt = define_casadi_symbolic(casadi_symbolic_mode,'N_vt',n_tangents);
                          case 'Lp'
                            p_vt = define_casadi_symbolic(casadi_symbolic_mode,'p_vt',n_tangents); % positive parts of tagnetial velocity (for switch detection)
                            n_vt = define_casadi_symbolic(casadi_symbolic_mode,'n_vt',n_tangents); % negative parts of tagnetial velocity (for switch detection)
                            alpha_vt = define_casadi_symbolic(casadi_symbolic_mode,'alpha_vt ',n_tangents); % step function of tangential velocities
                                                                                                            % impulse
                            P_vt = define_casadi_symbolic(casadi_symbolic_mode,'P_vt',n_tangents);
                            N_vt = define_casadi_symbolic(casadi_symbolic_mode,'N_vt',n_tangents);
                            Alpha_vt = define_casadi_symbolic(casadi_symbolic_mode,'Alpha_vt ',n_tangents);
                        end
                    end
                end
            end
            g_lift = [g_lift_theta_step;g_lift_beta];

            %% Collect algebaric varaibles for the specific DCS mode, define initial guess and bounds
            % TODO: @Anton: Do the bounds and guess specification already while defining the varaibles?

            z_impulse = []; % only for dcs_mode = cls, note that they are evaluated only at left boundary point of at FE
            switch dcs_mode
              case 'Stewart'
                % symbolic variables z = [theta;lambda;mu_Stewart];
                obj.z_all = [obj.theta;obj.lambda;obj.mu];
                z_switching = [obj.lambda;obj.mu];
              case 'Step'
                obj.z_all = [obj.alpha;obj.lambda_n;obj.lambda_p;obj.beta;obj.theta_step];
                z_switching = [obj.lambda_n;obj.lambda_p];
              case 'CLS'
                obj.z_all = [obj.lambda_normal;obj.y_gap];
                % Impulse
                %         z_impulse = [Lambda_normal;Y_gap;P_vn;N_vn];
                z_impulse = [obj.Lambda_normal;obj.Y_gap;L_vn];
                if friction_exists
                    % tangetial contact froce (firction force)
                    obj.z_all = [obj.z_all;obj.lambda_tangent];
                    % Impulse
                    z_impulse = [z_impulse;obj.Lambda_tangent]; % only for dcs_mode = cls, note that they are evaluated only at left boundary point of at FE
                    % friction aux multipliers
                    if isequal(friction_model,'Polyhedral')
                        % polyhedral friction model algebaric variables
                        obj.z_all = [obj.z_all;obj.gamma_d;obj.beta_d;obj.delta_d];
                        % Polyhedral friction - collect impulse variables
                        z_impulse = [z_impulse;obj.Gamma_d;obj.Beta_d;obj.Delta_d];
                    end
                    if isequal(friction_model,'Conic')
                        % conic friction model algebaric variables
                        obj.z_all = [obj.z_all;obj.gamma;obj.beta];
                        % Conic impulse
                        z_impulse = [z_impulse;obj.Gamma;obj.Beta];
                        switch conic_model_switch_handling
                          case 'Plain'
                            % no extra constraints
                          case 'Abs'
                            obj.z_all = [obj.z_all;obj.p_vt;obj.n_vt];
                            % Impulse
                            z_impulse = [z_impulse;obj.P_vt;obj.N_vt];
                          case 'Lp'
                            obj.z_all = [obj.z_all;obj.p_vt;obj.n_vt;obj.alpha_vt];
                            % Impulse
                            z_impulse = [z_impulse;obj.P_vt;obj.N_vt;obj.Alpha_vt];
                        end
                    end
                end
            end
            %% Add user provided algebraic
            obj.z_all = vertcat(obj.z_all,obj.z);
        end
        
        function verify_and_backfill(obj, settings)
            import casadi.*
            dims = obj.dims;
            if size(obj.x, 1) ~= 0
                dims.n_x = length(obj.x);
                % check  lbx
                if size(obj.lbx, 1) ~= 0
                    if length(obj.lbx) ~= dims.n_x
                        error('nosnoc: The vector lbx, for the lower bounds of x has the wrong size.')
                    end
                else
                    obj.lbx = -inf*ones(dims.n_x,1);
                end
                % check ubx
                if size(obj.ubx, 1) ~= 0
                    if length(obj.ubx) ~= dims.n_x
                        error('nosnoc: The vector ubx, for the upper bounds of x has the wrong size.')
                    end
                else
                    obj.ubx = inf*ones(dims.n_x,1);
                end
            else
                error('nosnoc: Please provide the state vector x, a CasADi symbolic variable.');
            end

            %% Check is u provided
            if size(obj.u, 1) ~= 0
                dims.n_u = length(obj.u);
                % check  lbu
                if size(obj.lbu, 1) ~= 0
                    if length(obj.lbu) ~= dims.n_u
                        error('nosnoc: The vector lbu, for the lower bounds of u has the wrong size.')
                    end
                else
                    obj.lbu = -inf*ones(dims.n_u,1);
                end
                % check ubu
                if size(obj.ubu, 1) ~= 0
                    if length(obj.ubu) ~= dims.n_u
                        error('nosnoc: The vector ubu, for the upper bounds of u has the wrong size.')
                    end
                else
                    obj.ubu = inf*ones(dims.n_u,1);
                end
                % check u0
                if size(obj.u0, 1) ~= 0
                    if length(obj.u0) ~= dims.n_u
                        error('nosnoc: The vector u0, for the initial guess of u has the wrong size.')
                    end
                else
                    obj.u0 = 0*ones(dims.n_u,1);
                end
            else
                obj.u = define_casadi_symbolic(settings.casadi_symbolic_mode,'',0);
                obj.u0 = [];
                dims.n_u = 0;
                if settings.print_level >=1
                    fprintf('nosnoc: No control vector u is provided. \n')
                end
                obj.lbu = [];
                obj.ubu = [];
            end
            %% Check if z is provided
            if size(obj.z, 1) ~= 0
                dims.n_z = length(obj.z);

                if size(obj.z0, 1) ~= 0
                    if length(obj.z0) ~= dims.n_z
                        error('nosnoc: The vector z0, for the initial guess of z has the wrong size.')
                    end
                else
                    obj.z0 = zeros(dims.n_z, 1);
                end

                if size(obj.lbz, 1) ~= 0
                    if length(obj.lbz) ~= dims.n_z
                        error('nosnoc: The vector lbz, for the lower bound of z has the wrong size.')
                    end
                else
                    obj.lbz = -inf*ones(dims.n_z, 1);
                end

                if size(obj.ubz, 1) ~= 0
                    if length(obj.ubz) ~= dims.n_z
                        error('nosnoc: The vector ubz, for the lower bound of z has the wrong size.')
                    end
                else
                    obj.ubz = inf*ones(dims.n_z, 1);
                end
            else
                dims.n_z = 0;
                obj.z0 = [];
                obj.lbz = [];
                obj.ubz = [];
                obj.z = define_casadi_symbolic(settings.casadi_symbolic_mode,'',0);
            end
            %% Global vars (i.e., variables that do not change with time)
            if size(obj.v_global, 1) ~= 0
                n_v_global = length(obj.v_global);
                if size(obj.v0_global, 1) ~= 0
                    if length(obj.v0_global) ~= n_v_global
                        error('nosnoc: The vector v0_global, for the initial guess of v_global has the wrong size.')
                    end
                else
                    obj.v0_global = zeros(n_v_global, 1);
                end

                if size(obj.lbv_global, 1) ~= 0
                    if length(obj.lbv_global) ~= n_v_global
                        error('nosnoc: The vector lbv_global, for the lower bound of v_global has the wrong size.')
                    end
                else
                    obj.lbv_global = -inf*ones(n_v_global, 1);
                end

                if size(obj.ubv_global, 1) ~= 0
                    if length(obj.ubv_global) ~= dims.n_v_global
                        error('nosnoc: The vector ubv_global, for the upper bound of v_global has the wrong size.')
                    end
                else
                    obj.ubv_global = inf*ones(n_v_global, 1);
                end
            else
                n_v_global = 0;
                obj.v_global = define_casadi_symbolic(settings.casadi_symbolic_mode, '', 0);
                obj.v0_global = [];
                obj.lbv_global = [];
                obj.ubv_global = [];
            end

            %% Parameters (time variable and that do not change with time)
            if size(obj.p_global, 1) ~= 0
                dims.n_p_global = size(obj.p_global,1);
                if size(obj.p_global_val, 1) ~= 0
                    if size(obj.p_global_val,1) ~= dims.n_p_global
                        error('nosnoc: User provided p_global_val has the wrong size.')
                    end
                else
                    p_global_val = zeros(dims.n_p_global,1);
                end
            else
                dims.n_p_global = 0;
                obj.p_global = define_casadi_symbolic(settings.casadi_symbolic_mode,'',0);
                obj.p_global_val = [];
                if settings.print_level >= 1
                    fprintf('nosnoc: No global parameters given. \n')
                end
            end

            if size(obj.p_time_var, 1) ~= 0
                dims.n_p_time_var = size(obj.p_time_var, 1);
                if size(obj.p_time_var_val, 1) ~= 0
                    if size(obj.p_time_var_val) ~= [dims.n_p_time_var, dims.N_stages]
                        error('nosnoc: User provided p_global_val has the wrong size.')
                    end
                else
                    obj.p_time_var_val = zeros(dims.n_p_time_var, dims.N_stages);
                end

                obj.p_time_var_stages = [];
                for ii=1:dims.N_stages
                    var_full = define_casadi_symbolic(settings.casadi_symbolic_mode, ['p_time_var_' num2str(ii)], dims.n_p_time_var);
                    obj.p_time_var_stages = horzcat(obj.p_time_var_stages, var_full);
                end
            else
                dims.n_p_time_var = 0;
                obj.p_time_var = define_casadi_symbolic(settings.casadi_symbolic_mode,'',0);
                obj.p_time_var_stages = define_casadi_symbolic(settings.casadi_symbolic_mode,'', [0, dims.N_stages]);
                obj.p_time_var_val = double.empty(0,dims.N_stages);
                if settings.print_level >= 1
                    fprintf('nosnoc: No time varying parameters given. \n')
                end
            end

            p = vertcat(obj.p_global,obj.p_time_var);

            %% Stage and terminal costs check
            if ~size(obj.f_q, 1) ~= 0
                if settings.print_level >=1
                    fprintf('nosnoc: No stage cost is provided. \n')
                end
                obj.f_q = 0;
            end

            if size(obj.f_q_T, 1) ~= 0
                terminal_cost = 1;
            else
                if settings.print_level >=1
                    fprintf('nosnoc: No terminal cost is provided. \n')
                end
                obj.f_q_T = 0;
            end
            %% Least squares objective terms with variables references
            if size(obj.lsq_x, 1) ~= 0
                if length(lsq_x)<3
                    error('nosnoc: In lsq_x either the least squares function, the reference of the weight matrix are missing.')
                end
                if size(lsq_x{2},1)~=size(lsq_x{1})
                    error('nosnoc: The dimensions of the least squares error term and weighting matrix for the differential states do not match.')
                end
                if size(lsq_x{1},1)~=size(lsq_x{3})
                    error('nosnoc: The dimensions of the least squares error term and reference for the differential states do not match.')
                end

                n_x_ref_rows = size(lsq_x{2},1);
                n_x_ref_cols = size(lsq_x{2},2);
                if n_x_ref_cols == dims.N_stages
                    fprintf('nosnoc: the provided reference for the differential states is time variable. \n');
                elseif n_x_ref_cols == 1
                    % replaciate
                    fprintf('nosnoc: the provided reference for the differential states is constant over time. \n');
                    lsq_x{2} = repmat(lsq_x{2},1,dims.N_stages);
                else
                    fprintf('nosnoc: The reference in lsq_x has to have a length of %d (if constant) or %d if time vriables. \n',1,dims.N_stages)
                    error('nosnoc: Please provide x_ref in lsq_x{1} with an appropaite size.')
                end
                x_ref_val = lsq_x{2};
                x_ref = define_casadi_symbolic(settings.casadi_symbolic_mode,'x_ref',n_x_ref_rows);
                f_lsq_x = (lsq_x{1}-x_ref)'*lsq_x{3}*(lsq_x{1}-x_ref);
            else
                x_ref = define_casadi_symbolic(settings.casadi_symbolic_mode,'x_ref',1);
                f_lsq_x = 0;
                x_ref_val = zeros(1,dims.N_stages);
            end

            % least square terms for control inputs
            if size(obj.lsq_u, 1) ~= 0
                if length(obj.lsq_u)<3
                    error('nosnoc: In lsq_u either the least squares function, the reference of the weight matrix are missing.')
                end
                if size(obj.lsq_u{2},1)~=size(obj.lsq_u{1})
                    error('nosnoc: The dimensions of the least squares error term and weighting matrix for the control input do not match.')
                end
                if size(obj.lsq_u{1},1)~=size(obj.lsq_u{3})
                    error('nosnoc: The dimensions of the least squares error term and reference for the control input do not match.')
                end
                n_u_ref_rows = size(obj.lsq_u{2},1);
                n_u_ref_cols = size(obj.lsq_u{2},2);
                if n_u_ref_cols == dims.N_stages
                    fprintf('nosnoc: the provided reference for the control inputs is time variable. \n');
                elseif n_u_ref_cols == 1
                    % replaciate
                    fprintf('nosnoc: the provided reference for the control inputs is constant over time. \n');
                    lsq_u{2} = repmat(lsq_u{2},1,dims.N_stages);
                else
                    fprintf('nosnoc: The reference in lsq_u has to have a length of %d (if constant) or %d if time vriables. \n',1,dims.N_stages)
                    error('nosnoc: Please provide u_ref in lsq_u{2} with an appropaite size.')
                end
                u_ref_val = lsq_u{2};
                u_ref = define_casadi_symbolic(settings.casadi_symbolic_mode,'u_ref',n_u_ref_rows);
                f_lsq_u = (lsq_u{1}-u_ref)'*lsq_u{3}*(lsq_u{1}-u_ref);
            else
                u_ref = define_casadi_symbolic(settings.casadi_symbolic_mode,'u_ref',1);
                f_lsq_u = 0;
                u_ref_val = zeros(1,dims.N_stages);
            end


            % least square terms for control inputs
            if size(obj.lsq_T, 1) ~= 0
                % sanity chkecs on the input
                if length(obj.lsq_T)<3
                    error('nosnoc: In lsq_u either the least squares function, the reference or the weight matrix are missing.')
                end
                if size(obj.lsq_T{2},1)~=size(obj.lsq_T{1})
                    error('nosnoc: The dimensions of the least squares error term and weighting matrix for the terminal cost do not match.')
                end
                if size(obj.lsq_T{1},1)~=size(obj.lsq_T{3})
                    error('nosnoc: The dimensions of the least squares error term and reference for the terminal cost do not match.')
                end

                n_x_T_rows = size(lsq_T{2},1);
                n_x_T_cols = size(lsq_T{2},2);
                if n_x_T_cols == 1
                    fprintf('nosnoc: the provided reference for the terminal cost is ok. \n');
                else
                    fprintf('nosnoc: The reference in lsq_T has to be a vector of length %d. \n',length(lsq_T{1}));
                    error('nosnoc: Please provide a reference vector in lsq_T{2} with an appropaite size.')
                end
                x_ref_end_val = lsq_T{2};
                x_ref_end = define_casadi_symbolic(settings.casadi_symbolic_mode,'x_ref_end',n_x_T_rows);
                f_lsq_T = (lsq_T{1}-x_ref_end)'*lsq_T{3}*(lsq_T{1}-x_ref_end);
            else
                x_ref_end  = define_casadi_symbolic(settings.casadi_symbolic_mode,'x_ref_end',1);
                f_lsq_T = 0;
                x_ref_end_val = 0;
            end

            %% Inequality constraints check
            if size(obj.g_path, 1) ~= 0
                g_path_constraint  = 1;
                n_g_path = length(g_path);
                if size(obj.g_path_lb, 1) ~= 0
                    if length(g_path_lb)~=n_g_path;
                        error('The user provided vector g_path_lb has the wrong size.')
                    end
                else
                    g_path_lb = -inf*ones(n_g_path,1);
                end

                if size(obj.g_path_ub, 1) ~= 0
                    if length(g_path_ub)~=n_g_path;
                        error('The user provided vector g_path_ub has the wrong size.')
                    end
                else
                    g_path_ub =  0*ones(n_g_path,1);
                end
                g_path_fun  = Function('g_path_fun',{x,u,p,v_global},{g_path});
            else
                n_g_path = 0;
                g_path_constraint  = 0;
                if settings.print_level >=1
                    fprintf('nosnoc: No path constraints are provided. \n')
                end
            end

            %% Check path complementarity constraints
            g_comp_path_constraint  = 0;
            if size(obj.g_comp_path, 1) ~= 0
                g_comp_path_constraint  = 1;
                if size(g_comp_path, 2) ~= 2
                    error('g_comp_path must be of size (m, 2)')
                end
                g_comp_path_fun  = Function('g_comp_path_fun',{x,u,p,v_global},{g_comp_path});
            else
                g_comp_path_constraint = 0;
                if settings.print_level >=1
                    fprintf('nosnoc: No path complementarity constraints are provided. \n')
                end
            end
            %% Terminal constraints
            if size(obj.g_terminal, 1) ~= 0
                terminal_constraint = 1;
                n_g_terminal = length(obj.g_terminal);
                if size(obj.g_terminal_lb, 1) ~= 0
                    if length(g_terminal_lb)~=n_g_terminal
                        error('nosnoc: The provided vector g_terminal_lb has the wrong size.')
                    end
                else
                    g_terminal_lb = 0*ones(n_g_terminal,1);
                end

                if size(obj.g_terminal_ub, 1) ~= 0
                    if length(g_terminal_ub)~=n_g_terminal
                        error('nosnoc: The provided vector g_terminal_ub has the wrong size.')
                    end
                else
                    g_terminal_ub =  0*ones(n_g_terminal,1);
                end
                g_terminal_fun  = Function('g_terminal_fun',{x,p_global,v_global},{g_terminal});
            else
                terminal_constraint = 0;
                n_g_terminal = 0;
                if settings.print_level >=1
                    fprintf('nosnoc: No terminal constraints are provided. \n')
                end
            end

            obj.g_Stewart = {};
            c_all = [];
            obj.friction_exists = 0;

            if isequal(settings.dcs_mode,'CLS')
                % TODO: there is some repetition to the time_freezing check, this should be unified!!!!
                % Check existence of relevant functions
                dims.n_sys = 1; % always one subystem in CLS (only loops over n_contacts later)
                if isempty(obj.f_c)
                    error('nosnoc: Please provide the gap functions model.f_c.')
                end
                dims.n_contacts = length(obj.f_c);

                % coefficient of friction checks
                if size(obj.mu, 1) ~= 0
                    if length(obj.mu) ~= 1 && length(obj.mu) ~= dims.n_contacts
                        error('The length of model.mu has to be one or match the length of model.f_c')
                    end
                    if length(obj.mu) == 1
                        mu = mu*ones(dims.n_contacts,1);
                        obj.mu = mu;
                    end

                    if any(obj.mu > 0)
                        obj.friction_exists = 1;
                    else
                        obj.friction_exists = 0;
                    end
                else
                    obj.mu = zeros(dims.n_contacts,1);
                    fprintf('nosnoc: Coefficients of friction mu not provided, setting it to zero for all contacts. \n')
                end
                if any(obj.mu<0)
                    error('nosnoc: The coefficients of friction mu should be nonnegative.')
                end

                % coefficent of restiution check
                if isempty(obj.e)
                    error('nosnoc:  Please provide a coefficient of restitution via model.e')
                else
                    if length(obj.e) ~= 1 && length(obj.e) ~= dims.n_contacts
                        error('The length of model.e has to be one or match the length of model.f_c')
                    end
                    if length(obj.e) == 1
                        e = e*ones(dims.n_contacts,1);
                        obj.e = e;
                    end
                end
                if any(abs(1-e)>1) || any(e<0)
                    error('nosnoc: the coefficient of restitution e should be in [0,1].')
                end

                % dimensions and state space split
                settings.casadi_symbolic_mode = obj.x(1).type_name();
                if mod(n_x,2)
                    dims.n_q = (n_x-1)/2;
                else
                    dims.n_q = n_x/2;
                end
                if isempty(obj.q) && isempty(obj.v)
                    q = x(1:dims.n_q);
                    v = x(dims.n_q+1:2*dims.n_q);
                end

                if isempty(obj.f_v)
                    error('nosnoc: the function f_v (collecting all generalized forces), in M(q) = dv/dt =  f_v(q,v,u) + J_n\lambda_n +J_t\lambda_t ~ is not provided in model.');
                end

                % Check intertia matrix
                if isempty(obj.M)
                    fprintf('nosnoc: Inertia matrix M is not provided, set to default value: M = eye(n_q). \n')
                    M = eye(dims.n_q);
                    invM = inv(M);
                else
                    invM = inv(M);
                end

                %  Normal Contact Jacobian
                if size(obj.J_normal, 1) ~= 0
                    J_normal = obj.J_normal;
                    J_normal_exists = 1;
                else
                    J_normal_exists = 0;
                end

                if J_normal_exists
                    if size(J_normal,1)~=dims.n_q && size(J_normal,2)~=dims.n_contacts
                        fprintf('nosnoc: J_normal should be %d x %d matrix.\n',dims.n_q,dims.n_contacts);
                        error('nosnoc: J_normal has the wrong size.')
                    end
                    J_normal_exists = 1;
                else
                    J_normal = f_c.jacobian(q)';
                    fprintf('nosnoc: normal contact Jacobian not provided, but it is computed from the gap functions.\n');
                    J_normal_exists = 1;
                end

                if is_zero(J_normal)
                    error('nosnoc: The normal vector should have at least one non-zero entry.')
                end

                % Tangent Contact Jacobian
                if obj.friction_exists
                    if isequal(settings.friction_model,'Conic')
                        if size(obj.J_tangent, 1) ~= 0
                            J_tangent = obj.J_tangent;
                            if size(J_tangent,1)~=dims.n_q
                                error('nosnoc: J_tangent has the wrong size.')
                            end
                        else
                            error('nosnoc: please provide the tangent Jacobian in model.J_tangent.')
                        end
                    end

                    if isequal(settings.friction_model,'Polyhedral')
                        if isempty(obj.D_tangent)
                            error('nosnoc: please provide the polyhedral tangent Jacobian in model.D_tangent, e.g., using the conic tangent Jacobian model.J_tangent: D_tangent = [J_tangent(q_0),-J_tangent(q_0)].')
                        end
                    end
                end
                % Dimension of tangents
                dims.n_t = 0;
                if obj.friction_exists
                    if isequal(friction_model,'Polyhedral')
                        dims.n_t = size(D_tangent,2)/dims.n_contacts; % number of tanget multipliers for a single contactl
                    elseif isequal(friction_model,'Conic')
                        dims.n_t = size(J_tangent,2)/dims.n_contacts; % number of tanget multipliers for a single contactl
                    end
                    dims.n_tangents = dims.n_t*dims.n_contacts; % number tangent forces for all multpliers
                else
                    dims.n_tangents = 0;
                end
            end

            if isequal(settings.dcs_mode,'Step') || isequal(settings.dcs_mode,'Stewart')
                if isempty(obj.F)
                    % Don't need F
                    if ~settings.general_inclusion
                        error('nosnoc: Matrix F (or matrices F_i) with PSS modes not provided.');
                    else
                        % TODO Implement more subsystems.
                        dims.n_sys = 1;
                    end
                else
                    % check how many subsystems are present
                    if iscell(obj.F)
                        dims.n_sys = length(obj.F);
                    else
                        obj.F = {obj.F};
                        dims.n_sys = 1;
                    end
                end

                if isempty(obj.S)
                    % if we are using general inclusions we dont need S.
                    if ~settings.general_inclusion
                        % if the matrix S is not provided, maybe the g_ind are available
                        % directly?
                        if isequal(settings.dcs_mode,'Stewart')
                            if exist('g_ind')
                                if ~iscell(g_ind)
                                    g_ind = {g_ind};
                                end

                                for ii = 1:dims.n_sys
                                    % discriminant functions
                                    obj.g_Stewart{ii} = g_ind{ii};
                                    c_all = [c_all; zeros(1,settings.casadi_symbolic_mode)];
                                end
                            else
                                error(['nosnoc: Neither the sign matrix S nor the indicator functions g_ind for regions are provided. ' ...
                                        'Either provide the matrix S and the expression for c, or the expression for g_ind.']);
                            end
                        else
                            error(['nosnoc: The user uses settings.dcs_mode = ''Step'', but the sign matrix S is not provided. Please provide the matrix S and the expressions for c(x) (definfing the region boundaries).']);
                        end
                    else
                        if isempty(obj.c)
                            error('nosnoc: Expresion for c, the constraint function for regions R_i is not provided.');
                        else
                            if ~iscell(obj,c)
                                obj.c = {obj.c};
                            end
                            if length(obj.c) ~= dims.n_sys
                                error('nosnoc: Number of different expressions for c does not match number of subsystems.')
                            end
                            for ii = 1:dims.n_sys
                                c_all = [c_all; c{ii}];
                                n_c{ii} = length(c{ii});
                                dims.n_c_sys  = [dims.n_c_sys;length(c{ii})];
                            end

                        end
                    end
                else
                    % Check if all data is avilable and if dimensions match.
                    if ~iscell(obj.S)
                        obj.S = {obj.S};
                    end
                    if length(obj.S) ~= dims.n_sys
                        error('nosnoc: Number of matrices S does not match number of subsystems. Note that the number of subsystems is taken to be number of matrices F_i which collect the modes of every subsystem.')
                    end
                    % Check constraint function c
                    if isempty(obj.c)
                        error('nosnoc: Expresion for c, the constraint function for regions R_i is not provided.');
                    else
                        if ~iscell(obj.c)
                            obj.c = {obj.c};
                        end
                        if length(obj.c) ~= dims.n_sys
                            error('nosnoc: Number of different expressions for c does not match number of subsystems (taken to be number of matrices F_i which collect the modes of every subsystem).')
                        end
                    end

                    % check are the matrices dense
                    if isequal(settings.dcs_mode,'Stewart')
                        for ii = 1:dims.n_sys
                            if any(sum(abs(obj.S{ii}),2)<size(obj.S{ii},2))
                                if dims.n_sys == 1
                                    error('nosnoc: The matrix S is not dense. Either provide a dense matrix or use settings.mode = ''Step''.');
                                else
                                    error(['The matrix S{' num2str(ii) '} of the provided matrices is not dense. Either provide all dense matrices or use settings.mode = ''Step''.']);
                                end
                            end
                        end
                    end

                    for ii = 1:dims.n_sys
                        if size(obj.S{ii},2) ~= length(obj.c{ii})
                            error('nosnoc: The matrix S and vector c do not have compatible dimension.');
                        end

                        % discrimnant functions
                        switch settings.dcs_mode
                            case 'Stewart'
                                % Create Stewart's indicator functions g_ind_ii
                                obj.g_Stewart{ii} = -obj.S{ii}*obj.c{ii};
                            case 'Step'
                                %eval(['c_' num2str(ii) '= c{ii};']);
                        end
                        % dimensions of c
                        c_all = [c_all; obj.c{ii}];
                        n_c{ii} = length(obj.c{ii});
                        dims.n_c_sys  = [dims.n_c_sys;length(obj.c{ii})];
                    end

                end

                if isempty(dims.n_c_sys)
                    dims.n_c_sys = 0;
                end

                if max(dims.n_c_sys) < 2 && isequal(settings.dcs_mode,'Step')
                    pss_lift_step_functions = 0;
                    if settings.print_level >=1
                        fprintf('nosnoc: settings.pss_lift_step_functions set to 0, as are step fucntion selections are already entering the ODE linearly.\n')
                    end
                end

                if ~settings.general_inclusion
                    dims.n_f_sys = arrayfun(@(sys) size(obj.F{sys},2),1:dims.n_sys);
                else
                    dims.n_f_sys = [size(f_x,1)];
                end
            end

            % populate functions that can already be generated
            obj.c_fun = Function('c_fun',{obj.x,p},{c_all});
            obj.g_Stewart_fun = Function('g_Stewart_fun',{obj.x,p},{vertcat(obj.g_Stewart{:})});
        end
    end  % methods
end % NosnocModel
