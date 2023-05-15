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
        lsq_x
        lsq_u
        lsq_T

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
            dims = NosnocDimensions();
        end

        function vars = generate_vars(obj,settings)

        end

        function reformulate(obj, settings)

        end

        function verify_and_backfill(obj, settings)
            if ~isempty(obj.x)
                obj.dims.n_x = length(x);
                % check  lbx
                if ~isempty(obj.lbx)
                    if length(model.lbx) ~= obj.dims.n_x
                        error('nosnoc: The vector lbx, for the lower bounds of x has the wrong size.')
                    end
                else
                    obj.lbx = -inf*ones(obj.dims.n_x,1);
                end
                % check ubx
                if ~isempty(obj.ubx)
                    if length(obj.ubx) ~= obj.dims.n_x
                        error('nosnoc: The vector ubx, for the upper bounds of x has the wrong size.')
                    end
                else
                    obj.ubx = inf*ones(obj.dims.n_x,1);
                end
            else
                error('nosnoc: Please provide the state vector x, a CasADi symbolic variable.');
            end

            %% Check is u provided
            if ~isempty(obj.u)
                obj.dims.n_u = length(obj.u);
                % check  lbu
                if ~isempty(obj.lbu)
                    if length(obj.lbu) ~= obj.dims.n_u
                        error('nosnoc: The vector lbu, for the lower bounds of u has the wrong size.')
                    end
                else
                    obj.lbu = -inf*ones(obj.dims.n_u,1);
                end
                % check ubu
                if ~isempty(obj.ubu)
                    if length(obj.ubu) ~= obj.dims.n_u
                        error('nosnoc: The vector ubu, for the upper bounds of u has the wrong size.')
                    end
                else
                    obj.ubu = inf*ones(obj.dims.n_u,1);
                end
                % check u0
                if ~isempty(obj.u0)
                    if length(u0) ~= obj.dims.n_u
                        error('nosnoc: The vector u0, for the initial guess of u has the wrong size.')
                    end
                else
                    obj.u0 = 0*ones(obj.dims.n_u,1);
                end
            else
                obj.u = define_casadi_symbolic(casadi_symbolic_mode,'',0);
                obj.u0 = [];
                obj.dims.n_u = 0;
                if print_level >=1
                    fprintf('nosnoc: No control vector u is provided. \n')
                end
                obj.lbu = [];
                obj.ubu = [];
            end
            %% Check if z is provided
            if ~isempty(obj.z)
                obj.dims.n_z = length(z);

                if ~isempty(obj.z0)
                    if length(obj.z0) ~= obj.dims.n_z
                        error('nosnoc: The vector z0, for the initial guess of z has the wrong size.')
                    end
                else
                    obj.z0 = zeros(obj.dims.n_z, 1);
                end

                if ~isempty(obj.lbz)
                    if length(obj.lbz) ~= obj.dims.n_z
                        error('nosnoc: The vector lbz, for the lower bound of z has the wrong size.')
                    end
                else
                    obj.lbz = -inf*ones(obj.dims.n_z, 1);
                end

                if ~isempty(obj.ubz)
                    if length(obj.ubz) ~= obj.dims.n_z
                        error('nosnoc: The vector ubz, for the lower bound of z has the wrong size.')
                    end
                else
                    obj.ubz = inf*ones(obj.dims.n_z, 1);
                end
            else
                obj.dims.n_z = 0;
                obj.z0 = [];
                obj.lbz = [];
                obj.ubz = [];
                obj.z = define_casadi_symbolic(casadi_symbolic_mode,'',0);
            end
            %% Global vars (i.e., variables that do not change with time)
            if ~isempty(obj.v_global)
                obj.obj.dims.n_v_global = length(obj.v_global);
                if ~isempty(obj.v0_global)
                    if length(obj.v0_global) ~= obj.dims.n_v_global
                        error('nosnoc: The vector v0_global, for the initial guess of v_global has the wrong size.')
                    end
                else
                    obj.v0_global = zeros(obj.dims.n_z, 1);
                end

                if ~isempty(obj.lbv_global)
                    if length(obj.lbv_global) ~= obj.dims.n_v_global
                        error('nosnoc: The vector lbv_global, for the lower bound of v_global has the wrong size.')
                    end
                else
                    obj.lbv_global = -inf*ones(obj.dims.n_z, 1);
                end

                if ~isempty(obj.ubv_global)
                    if length(obj.ubv_global) ~= obj.dims.n_v_global
                        error('nosnoc: The vector ubv_global, for the upper bound of v_global has the wrong size.')
                    end
                else
                    obj.ubv_global = inf*ones(obj.dims.n_z, 1);
                end
            else
                obj.dims.n_v_global = 0;
                obj.v_global = define_casadi_symbolic(casadi_symbolic_mode, '', 0);
                obj.v0_global = [];
                obj.lbv_global = [];
                obj.ubv_global = [];
            end

            %% Parameters (time variable and that do not change with time)
            if ~isempty(obj.p_global)
                obj.dims.n_p_global = size(obj.p_global,1);
                if ~isempty(obj.p_global_val)
                    if size(obj.p_global_val,1) ~= obj.dims.n_p_global
                        error('nosnoc: User provided p_global_val has the wrong size.')
                    end
                else
                    p_global_val = zeros(obj.dims.n_p_global,1);
                end
            else
                obj.dims.n_p_global = 0;
                obj.p_global = define_casadi_symbolic(casadi_symbolic_mode,'',0);
                obj.p_global_val = [];
                if print_level >= 1
                    fprintf('nosnoc: No global parameters given. \n')
                end
            end

            if ~isempty(obj.p_time_var)
                obj.dims.n_p_time_var = size(obj.p_time_var, 1);
                if ~isempty(obj.p_time_var_val)
                    if size(obj.p_time_var_val) ~= [obj.dims.n_p_time_var, N_stages]
                        error('nosnoc: User provided p_global_val has the wrong size.')
                    end
                else
                    obj.p_time_var_val = zeros(obj.dims.n_p_time_var, N_stages);
                end

                obj.p_time_var_stages = [];
                for ii=1:N_stages
                    var_full = define_casadi_symbolic(casadi_symbolic_mode, ['p_time_var_' num2str(ii)], obj.dims.n_p_time_var);
                    obj.p_time_var_stages = horzcat(obj.p_time_var_stages, var_full);
                end
            else
                obj.dims.n_p_time_var = 0;
                obj.p_time_var = define_casadi_symbolic(casadi_symbolic_mode,'',0);
                obj.p_time_var_stages = define_casadi_symbolic(casadi_symbolic_mode,'', [0, N_stages]);
                obj.p_time_var_val = double.empty(0,N_stages);
                if print_level >= 1
                    fprintf('nosnoc: No time varying parameters given. \n')
                end
            end

            p = vertcat(p_global,p_time_var);

            %% g_z: stage algebraic constraints
            % TODO  long term: split up model_reformulation to allow f_alg to use the rest of stage Z
            if ~isempty(obj.g_z)
                obj.dims.n_g_z = length(obj.g_z);
            else
                obj.g_z = [];
                obj.dims.n_g_z = 0;
            end
            %% Stage and terminal costs check
            if ~~isempty(obj.f_q)
                if print_level >=1
                    fprintf('nosnoc: No stage cost is provided. \n')
                end
                %     eval(['f_q = ', casadi_symbolic_mode, '.zeros(1);'])
                obj.f_q = 0;
            end

            if ~isempty(obj.f_q_T)
                terminal_cost = 1;
            else
                if print_level >=1
                    fprintf('nosnoc: No terminal cost is provided. \n')
                end
                %     eval(['f_q_T = ', casadi_symbolic_mode, '.zeros(1);'])
                obj.f_q_T = 0;
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
                n_g_terminal = length(obj.g_terminal);
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

            dcs_mode = settings.dcs_mode;
            g_Stewart = {};
            g_ind_vec = [];
            c_all = [];
            obj.dims.m_vec = [];
            n_c_sys = [];
            obj.dims.n_q = [];
            obj.friction_exists = 0;

            if isequal(dcs_mode,'CLS')
                % TODO: there is some repetition to the time_freezing check, this should be unified!!!!
                % Check existence of relevant functions
                obj.dims.n_sys = 1; % always one subystem in CLS (only loops over n_contacts later)
                if isempty(obj.f_c)
                    error('nosnoc: Please provide the gap functions model.f_c.')
                end
                obj.dims.n_contacts = length(obj.f_c);

                % coefficient of friction checks
                if ~isempty(obj.mu)
                    if length(obj.mu) ~= 1 && length(obj.mu) ~= obj.dims.n_contacts
                        error('The length of model.mu has to be one or match the length of model.f_c')
                    end
                    if length(obj.mu) == 1
                        mu = mu*ones(obj.dims.n_contacts,1);
                        obj.mu = mu;
                    end

                    if any(obj.mu > 0)
                        obj.friction_exists = 1;
                    else
                        obj.friction_exists = 0;
                    end
                else
                    obj.mu = zeros(obj.dims.n_contacts,1);
                    fprintf('nosnoc: Coefficients of friction mu not provided, setting it to zero for all contacts. \n')
                end
                if any(obj.mu<0)
                    error('nosnoc: The coefficients of friction mu should be nonnegative.')
                end

                % coefficent of restiution check
                if isempty(obj.e)
                    error('nosnoc:  Please provide a coefficient of restitution via model.e')
                else
                    if length(obj.e) ~= 1 && length(obj.e) ~= obj.dims.n_contacts
                        error('The length of model.e has to be one or match the length of model.f_c')
                    end
                    if length(obj.e) == 1
                        e = e*ones(obj.dims.n_contacts,1);
                        obj.e = e;
                    end
                end
                if any(abs(1-e)>1) || any(e<0)
                    error('nosnoc: the coefficient of restitution e should be in [0,1].')
                end

                % dimensions and state space split
                casadi_symbolic_mode = obj.x(1).type_name();
                if mod(n_x,2)
                    obj.dims.n_q = (n_x-1)/2;
                else
                    obj.dims.n_q = n_x/2;
                end
                if isempty(obj.q) && isempty(obj.v)
                    q = x(1:obj.dims.n_q);
                    v = x(obj.dims.n_q+1:2*obj.dims.n_q);
                end

                if isempty(obj.f_v)
                    error('nosnoc: the function f_v (collecting all generalized forces), in M(q) = dv/dt =  f_v(q,v,u) + J_n\lambda_n +J_t\lambda_t ~ is not provided in model.');
                end

                % Check intertia matrix
                if isempty(obj.M)
                    fprintf('nosnoc: Inertia matrix M is not provided, set to default value: M = eye(n_q). \n')
                    M = eye(obj.dims.n_q);
                    invM = inv(M);
                else
                    invM = inv(M);
                end

                %  Normal Contact Jacobian
                if ~isempty(obj.J_normal)
                    J_normal = obj.J_normal;
                    J_normal_exists = 1;
                else
                    J_normal_exists = 0;
                end

                if J_normal_exists
                    if size(J_normal,1)~=obj.dims.n_q && size(J_normal,2)~=obj.dims.n_contacts
                        fprintf('nosnoc: J_normal should be %d x %d matrix.\n',obj.dims.n_q,obj.dims.n_contacts);
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
                        if ~isempty(obj.J_tangent)
                            J_tangent = obj.J_tangent;
                            if size(J_tangent,1)~=obj.dims.n_q
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
                obj.dims.n_t = 0;
                if obj.friction_exists
                    if isequal(friction_model,'Polyhedral')
                        obj.dims.n_t = size(D_tangent,2)/obj.dims.n_contacts; % number of tanget multipliers for a single contactl
                    elseif isequal(friction_model,'Conic')
                        obj.dims.n_t = size(J_tangent,2)/obj.dims.n_contacts; % number of tanget multipliers for a single contactl
                    end
                    obj.dims.n_tangents = obj.dims.n_t*obj.dims.n_contacts; % number tangent forces for all multpliers
                else
                    obj.dims.n_tangents = 0;
                end
            end

            if isequal(dcs_mode,'Step') || isequal(dcs_mode,'Stewart')
                if isempty(obj.F)
                    % Don't need F
                    if ~settings.general_inclusion
                        error('nosnoc: Matrix F (or matrices F_i) with PSS modes not provided.');
                    else
                        % TODO Implement more subsystems.
                        obj.dims.n_sys = 1;
                        obj.dims.m_vec = [size(f_x,1)];
                    end
                else
                    % check how many subsystems are present
                    if iscell(obj.F)
                        obj.dims.n_sys = length(obj.F);
                    else
                        obj.F = {obj.F};
                        obj.dims.n_sys = 1;
                    end
                    % extract dimensions of subystems
                    for ii = 1:obj.dims.n_sys
                        m_temp = size(obj.F{ii},2);
                        obj.dims.m_vec  = [obj.dims.m_vec m_temp];
                    end
                end

                if isempty(obj.S)
                    % if we are using general inclusions we dont need S.
                    if ~settings.general_inclusion
                        % if the matrix S is not provided, maybe the g_ind are available
                        % directly?
                        if isequal(dcs_mode,'Stewart')
                            if exist('g_ind')
                                if ~iscell(g_ind)
                                    g_ind = {g_ind};
                                end

                                for ii = 1:obj.dims.n_sys
                                    % discriminant functions
                                    g_ind_vec =  [g_ind_vec;g_ind{ii};];
                                    g_Stewart{ii} = g_ind{ii};
                                    c_all = [c_all; zeros(1,casadi_symbolic_mode)];
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
                            if ~iscell(c)
                                c = {c};
                            end
                            if length(c) ~= obj.dims.n_sys
                                error('nosnoc: Number of different expressions for c does not match number of subsystems.')
                            end
                            for ii = 1:obj.dims.n_sys
                                c_all = [c_all; c{ii}];
                                n_c{ii} = length(c{ii});
                                n_c_sys  = [n_c_sys;length(c{ii})];
                            end

                        end
                    end
                else
                    % Check if all data is avilable and if dimensions match.
                    if ~iscell(obj.S)
                        obj.S = {obj.S};
                    end
                    if length(obj.S) ~= obj.dims.n_sys
                        error('nosnoc: Number of matrices S does not match number of subsystems. Note that the number of subsystems is taken to be number of matrices F_i which collect the modes of every subsystem.')
                    end
                    % Check constraint function c
                    if isempty(obj.c)
                        error('nosnoc: Expresion for c, the constraint function for regions R_i is not provided.');
                    else
                        if ~iscell(obj.c)
                            obj.c = {obj.c};
                        end
                        if length(c) ~= obj.dims.n_sys
                            error('nosnoc: Number of different expressions for c does not match number of subsystems (taken to be number of matrices F_i which collect the modes of every subsystem).')
                        end
                    end

                    % check are the matrices dense
                    if isequal(dcs_mode,'Stewart')
                        for ii = 1:obj.dims.n_sys
                            if any(sum(abs(obj.S{ii}),2)<size(obj.S{ii},2))
                                if obj.dims.n_sys == 1
                                    error('nosnoc: The matrix S is not dense. Either provide a dense matrix or use settings.mode = ''Step''.');
                                else
                                    error(['The matrix S{' num2str(ii) '} of the provided matrices is not dense. Either provide all dense matrices or use settings.mode = ''Step''.']);
                                end
                            end
                        end
                    end

                    for ii = 1:obj.dims.n_sys
                        if size(obj.S{ii},2) ~= length(c{ii})
                            error('nosnoc: The matrix S and vector c do not have compatible dimension.');
                        end

                        % discrimnant functions
                        switch dcs_mode
                            case 'Stewart'
                                % Create Stewart's indicator functions g_ind_ii
                                g_Stewart{ii} = -obj.S{ii}*c{ii};
                                g_ind_vec = [g_ind_vec ;-obj.S{ii}*c{ii}];
                            case 'Step'
                                %eval(['c_' num2str(ii) '= c{ii};']);
                        end
                        % dimensions of c
                        c_all = [c_all; c{ii}];
                        n_c{ii} = length(c{ii});
                        n_c_sys  = [n_c_sys;length(c{ii})];
                    end

                end
                % index sets and dimensions for ubsystems
                m_ind_vec = 1;
                if isequal(dcs_mode,'Step')
                    % double the size of the vectors, since alpha, 1-alpha treated at same time;
                    obj.dims.m_vec = sum(n_c_sys)*2;
                end
                for ii = 1:length(obj.dims.m_vec)-1
                    m_ind_vec  = [m_ind_vec,m_ind_vec(end)+obj.dims.m_vec(ii)];
                end
                % m_ind_vec = [cumsum(obj.dims.m_vec)-obj.dims.m_vec(1)+1]; % index ranges of the corresponding thetas and lambdas
                m = sum(obj.dims.m_vec);

                if isempty(n_c_sys)
                    n_c_sys = 0;
                end

                if max(n_c_sys) < 2 && isequal(dcs_mode,'Step')
                    pss_lift_step_functions = 0;
                    if print_level >=1
                        fprintf('nosnoc: settings.pss_lift_step_functions set to 0, as are step fucntion selections are already entering the ODE linearly.\n')
                    end
                end

                if ~settings.general_inclusion
                    obj.dims.n_f_sys = arrayfun(@(sys) size(obj.F{sys},2),1:obj.dims.n_sys);
                else
                    obj.dims.n_f_sys = [size(f_x,1)];
                end
            end
        end
        
    end

end
