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
            dims = NosnocDimensions();
        end

        function vars = generate_vars(obj,settings)

        end

        function reformulate(obj, settings)

        end

        function verify_and_backfill(obj, settings)
            if ~isempty(obj.x)
                n_x = length(x);
                % check  lbx
                if ~isempty(obj.lbx)
                    if length(model.lbx) ~= n_x
                        error('nosnoc: The vector lbx, for the lower bounds of x has the wrong size.')
                    end
                else
                    lbx = -inf*ones(n_x,1);
                end
                % check ubx
                if ~isempty(obj.ubx)
                    if length(model.ubx) ~= n_x
                        error('nosnoc: The vector ubx, for the upper bounds of x has the wrong size.')
                    end
                else
                    ubx = inf*ones(n_x,1);
                end
            else
                error('nosnoc: Please provide the state vector x, a CasADi symbolic variable.');
            end

            if ~isempty(obj.x)
                
            end
            %% Check is u provided
            if ~isempty(obj.u)
                n_u = length(model.u);
                % check  lbu
                if ~isempty(obj.lbu)
                    if length(model.lbu) ~= n_u
                        error('nosnoc: The vector lbu, for the lower bounds of u has the wrong size.')
                    end
                else
                    lbu = -inf*ones(n_u,1);
                end
                % check ubu
                if ~isempty(obj.ubu)
                    if length(model.ubu) ~= n_u
                        error('nosnoc: The vector ubu, for the upper bounds of u has the wrong size.')
                    end
                else
                    ubu = inf*ones(n_u,1);
                end
                % check u0
                if ~isempty(obj.u0)
                    if length(u0) ~= n_u
                        error('nosnoc: The vector u0, for the initial guess of u has the wrong size.')
                    end
                else
                    u0 = 0*ones(n_u,1);
                end
            else
                u = define_casadi_symbolic(casadi_symbolic_mode,'',0);
                u0 = [];
                n_u = 0;
                if print_level >=1
                    fprintf('nosnoc: No control vector u is provided. \n')
                end
                lbu = [];
                ubu = [];
            end
            %% Check if z is provided
            if ~isempty(obj.z)
                n_z = length(z);

                if ~isempty(obj.z0)
                    if length(model.z0) ~= n_z
                        error('nosnoc: The vector z0, for the initial guess of z has the wrong size.')
                    end
                else
                    z0 = zeros(n_z, 1);
                end

                if ~isempty(obj.lbz)
                    if length(model.lbz) ~= n_z
                        error('nosnoc: The vector lbz, for the lower bound of z has the wrong size.')
                    end
                else
                    lbz = -inf*ones(n_z, 1);
                end

                if ~isempty(obj.ubz)
                    if length(model.ubz) ~= n_z
                        error('nosnoc: The vector ubz, for the lower bound of z has the wrong size.')
                    end
                else
                    ubz = inf*ones(n_z, 1);
                end
            else
                n_z = 0;
                z0 = [];
                lbz = [];
                ubz = [];
                z = define_casadi_symbolic(casadi_symbolic_mode,'',0);
            end
            %% Global vars (i.e., variables that do not change with time)
            if ~isempty(obj.v_global)
                n_v_global = length(model.v_global);
                if ~isempty(obj.v0_global)
                    if length(model.v0_global) ~= n_v_global
                        error('nosnoc: The vector v0_global, for the initial guess of v_global has the wrong size.')
                    end
                else
                    z0 = zeros(n_z, 1);
                end

                if ~isempty(obj.lbv_global)
                    if length(model.lbv_global) ~= n_v_global
                        error('nosnoc: The vector lbv_global, for the lower bound of v_global has the wrong size.')
                    end
                else
                    lbz = -inf*ones(n_z, 1);
                end

                if ~isempty(obj.ubv_global)
                    if length(model.ubv_global) ~= n_v_global
                        error('nosnoc: The vector ubv_global, for the upper bound of v_global has the wrong size.')
                    end
                else
                    ubz = inf*ones(n_z, 1);
                end
            else
                n_v_global = 0;
                v_global = define_casadi_symbolic(casadi_symbolic_mode, '', 0);
                v0_global = [];
                lbv_global = [];
                ubv_global = [];
            end

            %% Parameters (time variable and that do not change with time)
            if ~isempty(obj.p_global)
                n_p_global = size(model.p_global,1);
                if ~isempty(obj.p_global_val)
                    if size(model.p_global_val,1) ~= n_p_global
                        error('nosnoc: User provided p_global_val has the wrong size.')
                    end
                else
                    p_global_val = zeros(n_p_global,1);
                end
            else
                n_p_global = 0;
                p_global = define_casadi_symbolic(casadi_symbolic_mode,'',0);
                p_global_val = [];
                if print_level >= 1
                    fprintf('nosnoc: No global parameters given. \n')
                end
            end

            if ~isempty(obj.p_time_var)
                n_p_time_var = size(model.p_time_var, 1);
                if ~isempty(obj.p_time_var_val)
                    if size(model.p_time_var_val) ~= [n_p_time_var, N_stages]
                        error('nosnoc: User provided p_global_val has the wrong size.')
                    end
                else
                    p_time_var_val = zeros(n_p_time_var, N_stages);
                end

                p_time_var_stages = [];
                for ii=1:N_stages
                    var_full = define_casadi_symbolic(casadi_symbolic_mode, ['p_time_var_' num2str(ii)], n_p_time_var);
                    p_time_var_stages = horzcat(p_time_var_stages, var_full);
                end
            else
                n_p_time_var = 0;
                p_time_var = define_casadi_symbolic(casadi_symbolic_mode,'',0);
                p_time_var_stages = define_casadi_symbolic(casadi_symbolic_mode,'', [0, N_stages]);
                p_time_var_val = double.empty(0,N_stages);
                if print_level >= 1
                    fprintf('nosnoc: No time varying parameters given. \n')
                end
            end

            p = vertcat(p_global,p_time_var);

            %% g_z: stage algebraic constraints
            % TODO  long term: split up model_reformulation to allow f_alg to use the rest of stage Z
            if ~isempty(obj.g_z)
                n_g_z = length(model.g_z);
            else
                g_z = [];
                n_g_z = 0;
            end
            %% Stage and terminal costs check
            if ~~isempty(obj.f_q)
                if print_level >=1
                    fprintf('nosnoc: No stage cost is provided. \n')
                end
                %     eval(['f_q = ', casadi_symbolic_mode, '.zeros(1);'])
                f_q = 0;
            end

            if ~isempty(obj.f_q_T)
                terminal_cost = 1;
            else
                if print_level >=1
                    fprintf('nosnoc: No terminal cost is provided. \n')
                end
                %     eval(['f_q_T = ', casadi_symbolic_mode, '.zeros(1);'])
                f_q_T = 0;
            end
            %% Least squares objective terms with variables references
            if ~isempty(obj.lsq_x)
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
                if n_x_ref_cols == N_stages
                    fprintf('nosnoc: the provided reference for the differential states is time variable. \n');
                elseif n_x_ref_cols == 1
                    % replaciate
                    fprintf('nosnoc: the provided reference for the differential states is constant over time. \n');
                    lsq_x{2} = repmat(lsq_x{2},1,N_stages);
                else
                    fprintf('nosnoc: The reference in lsq_x has to have a length of %d (if constant) or %d if time vriables. \n',1,N_stages)
                    error('nosnoc: Please provide x_ref in lsq_x{1} with an appropaite size.')
                end
                x_ref_val = lsq_x{2};
                x_ref = define_casadi_symbolic(casadi_symbolic_mode,'x_ref',n_x_ref_rows);
                f_lsq_x = (lsq_x{1}-x_ref)'*lsq_x{3}*(lsq_x{1}-x_ref);
            else
                x_ref = define_casadi_symbolic(casadi_symbolic_mode,'x_ref',1);
                f_lsq_x = 0;
                x_ref_val = zeros(1,N_stages);
            end

            % least square terms for control inputs
            if ~isempty(obj.lsq_u)
                if length(model.lsq_u)<3
                    error('nosnoc: In lsq_u either the least squares function, the reference of the weight matrix are missing.')
                end
                if size(model.lsq_u{2},1)~=size(model.lsq_u{1})
                    error('nosnoc: The dimensions of the least squares error term and weighting matrix for the control input do not match.')
                end
                if size(model.lsq_u{1},1)~=size(model.lsq_u{3})
                    error('nosnoc: The dimensions of the least squares error term and reference for the control input do not match.')
                end
                n_u_ref_rows = size(model.lsq_u{2},1);
                n_u_ref_cols = size(model.lsq_u{2},2);
                if n_u_ref_cols == N_stages
                    fprintf('nosnoc: the provided reference for the control inputs is time variable. \n');
                elseif n_u_ref_cols == 1
                    % replaciate
                    fprintf('nosnoc: the provided reference for the control inputs is constant over time. \n');
                    lsq_u{2} = repmat(lsq_u{2},1,N_stages);
                else
                    fprintf('nosnoc: The reference in lsq_u has to have a length of %d (if constant) or %d if time vriables. \n',1,N_stages)
                    error('nosnoc: Please provide u_ref in lsq_u{2} with an appropaite size.')
                end
                u_ref_val = lsq_u{2};
                u_ref = define_casadi_symbolic(casadi_symbolic_mode,'u_ref',n_u_ref_rows);
                f_lsq_u = (lsq_u{1}-u_ref)'*lsq_u{3}*(lsq_u{1}-u_ref);
            else
                u_ref = define_casadi_symbolic(casadi_symbolic_mode,'u_ref',1);
                f_lsq_u = 0;
                u_ref_val = zeros(1,N_stages);
            end


            % least square terms for control inputs
            if ~isempty(obj.lsq_T)
                % sanity chkecs on the input
                if length(model.lsq_T)<3
                    error('nosnoc: In lsq_u either the least squares function, the reference or the weight matrix are missing.')
                end
                if size(model.lsq_T{2},1)~=size(model.lsq_T{1})
                    error('nosnoc: The dimensions of the least squares error term and weighting matrix for the terminal cost do not match.')
                end
                if size(model.lsq_T{1},1)~=size(model.lsq_T{3})
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
                x_ref_end = define_casadi_symbolic(casadi_symbolic_mode,'x_ref_end',n_x_T_rows);
                f_lsq_T = (lsq_T{1}-x_ref_end)'*lsq_T{3}*(lsq_T{1}-x_ref_end);
            else
                x_ref_end  = define_casadi_symbolic(casadi_symbolic_mode,'x_ref_end',1);
                f_lsq_T = 0;
                x_ref_end_val = 0;
            end

            %% Inequality constraints check
            if ~isempty(obj.g_path)
                g_path_constraint  = 1;
                n_g_path = length(g_path);
                if ~isempty(obj.g_path_lb)
                    if length(g_path_lb)~=n_g_path;
                        error('The user provided vector g_path_lb has the wrong size.')
                    end
                else
                    g_path_lb = -inf*ones(n_g_path,1);
                end

                if ~isempty(obj.g_path_ub)
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
                if print_level >=1
                    fprintf('nosnoc: No path constraints are provided. \n')
                end
            end

            %% Check path complementarity constraints
            g_comp_path_constraint  = 0;
            if ~isempty(obj.g_comp_path)
                g_comp_path_constraint  = 1;
                if size(g_comp_path, 2) ~= 2
                    error('g_comp_path must be of size (m, 2)')
                end
                g_comp_path_fun  = Function('g_comp_path_fun',{x,u,p,v_global},{g_comp_path});
            else
                g_comp_path_constraint = 0;
                if print_level >=1
                    fprintf('nosnoc: No path complementarity constraints are provided. \n')
                end
            end
            %% Terminal constraints
            if ~isempty(obj.g_terminal)
                terminal_constraint = 1;
                n_g_terminal = length(model.g_terminal);
                if ~isempty(obj.g_terminal_lb)
                    if length(g_terminal_lb)~=n_g_terminal
                        error('nosnoc: The provided vector g_terminal_lb has the wrong size.')
                    end
                else
                    g_terminal_lb = 0*ones(n_g_terminal,1);
                end

                if ~isempty(obj.g_terminal_ub)
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
                if print_level >=1
                    fprintf('nosnoc: No terminal constraints are provided. \n')
                end
            end
        end
        
    end

end
