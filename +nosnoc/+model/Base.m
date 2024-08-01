classdef Base < matlab.mixin.Scalar & handle & matlab.mixin.CustomDisplay
% Base class for all ``nosnoc`` models. It contains shared properties such as the state, user algebraics, controls etc.
% It also contains the fields used to populate an Optimal Control problem such as Lagrange and Mayer cost terms,
% an interface for least squares costs as well as non-box path and terminal constraints.
    properties
        x % casadi.SX|casadi.MX: Differential state $x \in \mathbb{R}^{n_x}$
        lbx % double: $\underline{x} \in \mathbb{R}^{n_x}$, differential state lower bound.
        ubx % double: $\bar{x} \in \mathbb{R}^{n_x}$, differential state upper bound.
        x0 % double: $x_0 \in \mathbb{R}^{n_x}$, initial differential state, also used to initialize all differential state variables in the resulting MPCC.

        z % casadi.SX|casadi.MX: User algebraics $z \in \mathbb{R}^{n_z}$
        z0 % double: $z_0 \in \mathbb{R}^{n_z}$, used to initialize all user algebraic variables in the resulting MPCC.
        lbz % double: $\underline{z} \in \mathbb{R}^{n_z}$, user algebraic lower bound.
        ubz % double: $\bar{z} \in \mathbb{R}^{n_z}$, user algebraic upper bound.
        g_z % casadi.SX|casadi.MX: Constraint expression used to define the behavior of user algebraics.

        u % casadi.SX|casadi.MX: Controls $u \in \mathbb{R}^{n_u}$.
        lbu % double: $\underline{u} \in \mathbb{R}^{n_u}$, controls lower bound.
        ubu % double: $\bar{u} \in \mathbb{R}^{n_u}$, controls upper bound.
        u0 % double: $u_0 \in \mathbb{R}^{n_u}$, used to initialize all control variables in the resulting MPCC.

        v_global % casadi.SX|casadi.MX: $\nu \in \mathbb{R}^{n_{\nu}}$ global variables (not time dependent).
        v0_global % double: $\nu_0 \in \mathbb{R}^{n_{\nu}}$, used to initialize all global variables in the resulting MPCC.
        lbv_global % double: $\underline{\nu} \in \mathbb{R}^{n_{\nu}}$, global variables lower bound.
        ubv_global % double: $\bar{\nu} \in \mathbb{R}^{n_{\nu}}$, global variables upper bound.
        
        p_global % casadi.SX|casadi.MX: Global parameters.
        p_global_val % double: Values for global parameters

        p_time_var % casadi.SX|casadi.MX: Time varying parameters which are considered to be constant over each control/integration interval.
        p_time_var_val % double: Values for time varying parameters.
        p % casadi.SX|casadi.MX: All model parameters

        f_q % casadi.SX|casadi.MX: Lagrange term cost.
        f_q_T % casadi.SX|casadi.MX: Mayer term cost.

        % least squares (TODO @Anton: more precise description, e.g. where are the weight matrices)
        % TODO(@anton) perhaps the calculated parts of the lsq interface should be read only

        lsq_x % cell: TODO describe
        x_ref % casadi.SX|casadi.MX:
        f_lsq_x % casadi.SX|casadi.MX:
        x_ref_val % double: 
        lsq_u % casadi.SX|casadi.MX:
        u_ref % casadi.SX|casadi.MX:
        f_lsq_u % casadi.SX|casadi.MX:
        u_ref_val % double: vector
        lsq_T % casadi.SX|casadi.MX:
        x_ref_end % casadi.SX|casadi.MX:
        f_lsq_T % casadi.SX|casadi.MX:
        x_ref_end_val % double: vector

        g_path % casadi.SX|casadi.MX: Path constraints.
        lbg_path % double: Lower bound on path constraints.
        ubg_path % double: Upper bound on path constraints.

        g_terminal % casadi.SX|casadi.MX: Terminal constraints.
        lbg_terminal % double: Lower bound on path constraints.
        ubg_terminal % double: Upper bound on path constraints.

        G_path % casadi.SX|casadi.MX: One half of path complementarities.
        H_path % casadi.SX|casadi.MX: One half of path complementarities.

        dims % struct: Dimensions struct, the contents of which depends on the subclass.
    end

    methods
        function obj = Base()
        end
        
        function verify_and_backfill(obj, opts)
            arguments
                obj
                opts nosnoc.Options
            end
            import casadi.*
            dims = obj.dims;

            if size(obj.x, 1) ~= 0
                dims.n_x = length(obj.x);
                % check  lbx
                if size(obj.lbx, 1) ~= 0
                    if length(obj.lbx) ~= dims.n_x
                        error('nosnoc:model:Base:lbx_size',...
                            'nosnoc: The vector lbx, for the lower bounds of x has the wrong size.')
                    end
                else
                    obj.lbx = -inf*ones(dims.n_x,1);
                end
                % check ubx
                if size(obj.ubx, 1) ~= 0
                    if length(obj.ubx) ~= dims.n_x
                        error('nosnoc:model:Base:ubx_size',...
                            'nosnoc: The vector ubx, for the upper bounds of x has the wrong size.')
                    end
                else
                    obj.ubx = inf*ones(dims.n_x,1);
                end
            else
                error('nosnoc:model:Base:x_missing',...
                    'nosnoc: Please provide the state vector x, a CasADi symbolic variable.');
            end
            opts.casadi_symbolic_mode = ['casadi.' obj.x(1).type_name()]; % TODO(@anton) this shoud probably live in opts class

            %% Check is u provided
            if size(obj.u, 1) ~= 0
                dims.n_u = length(obj.u);
                % check  lbu
                if size(obj.lbu, 1) ~= 0
                    if length(obj.lbu) ~= dims.n_u
                        error('nosnoc:model:Base:lbu_size',...
                            'nosnoc: The vector lbu, for the lower bounds of u has the wrong size.')
                    end
                else
                    obj.lbu = -inf*ones(dims.n_u,1);
                end
                % check ubu
                if size(obj.ubu, 1) ~= 0
                    if length(obj.ubu) ~= dims.n_u
                        error('nosnoc:model:Base:ubu_size',...
                            'nosnoc: The vector ubu, for the upper bounds of u has the wrong size.')
                    end
                else
                    obj.ubu = inf*ones(dims.n_u,1);
                end
                % check u0
                if size(obj.u0, 1) ~= 0
                    if length(obj.u0) ~= dims.n_u
                        error('nosnoc:model:Base:u0_size',...
                            'nosnoc: The vector u0, for the initial guess of u has the wrong size.')
                    end
                else
                    obj.u0 = 0*ones(dims.n_u,1);
                end
            else
                obj.u = define_casadi_symbolic(opts.casadi_symbolic_mode,'',0);
                dims.n_u = 0;
                obj.u0 = [];
                obj.lbu = [];
                obj.ubu = [];
            end
            
            %% Check if z is provided
            if size(obj.z, 1) ~= 0
                dims.n_z = length(obj.z);
                if size(obj.z0, 1) ~= 0
                    if length(obj.z0) ~= dims.n_z
                        error('nosnoc:model:Base:z0_size',...
                            'nosnoc: The vector z0, for the initial guess of z has the wrong size.')
                    end
                else
                    obj.z0 = zeros(dims.n_z, 1);
                end

                if size(obj.lbz, 1) ~= 0
                    if length(obj.lbz) ~= dims.n_z
                        error('nosnoc:model:Base:lbz_size',...
                            'nosnoc: The vector lbz, for the lower bound of z has the wrong size.')
                    end
                else
                    obj.lbz = -inf*ones(dims.n_z, 1);
                end

                if size(obj.ubz, 1) ~= 0
                    if length(obj.ubz) ~= dims.n_z
                        error('nosnoc:model:Base:ubz_size',...
                            'nosnoc: The vector ubz, for the lower bound of z has the wrong size.')
                    end
                else
                    obj.ubz = inf*ones(dims.n_z, 1);
                end
            else
                obj.z = define_casadi_symbolic(opts.casadi_symbolic_mode,'',0);
                dims.n_z = 0;
                obj.z0 = [];
                obj.lbz = [];
                obj.ubz = [];
            end
            
            %% Global vars (i.e., variables that do not change with time)
            if size(obj.v_global, 1) ~= 0
                dims.n_v_global = length(obj.v_global);
                if size(obj.v0_global, 1) ~= 0
                    if length(obj.v0_global) ~= dims.n_v_global
                        error('nosnoc:model:Base:v0_global_size',...
                            'nosnoc: The vector v0_global, for the initial guess of v_global has the wrong size.')
                    end
                else
                    obj.v0_global = zeros(dims.n_v_global, 1);
                end

                if size(obj.lbv_global, 1) ~= 0
                    if length(obj.lbv_global) ~= dims.n_v_global
                        error('nosnoc:model:Base:lbv_global_size',...
                            'nosnoc: The vector lbv_global, for the lower bound of v_global has the wrong size.')
                    end
                else
                    obj.lbv_global = -inf*ones(dims.n_v_global, 1);
                end

                if size(obj.ubv_global, 1) ~= 0
                    if length(obj.ubv_global) ~= dims.n_v_global
                        error('nosnoc:model:Base:ubv_global_size',...
                            'nosnoc: The vector ubv_global, for the upper bound of v_global has the wrong size.')
                    end
                else
                    obj.ubv_global = inf*ones(dims.n_v_global, 1);
                end
            else
                obj.v_global = define_casadi_symbolic(opts.casadi_symbolic_mode, '', 0);
                dims.n_v_global = 0;
                obj.v0_global = [];
                obj.lbv_global = [];
                obj.ubv_global = [];
            end

            %% Parameters (time variable and that do not change with time)
            if size(obj.p_global, 1) ~= 0
                dims.n_p_global = size(obj.p_global,1);
                if size(obj.p_global_val, 1) ~= 0
                    if size(obj.p_global_val,1) ~= dims.n_p_global
                        error('nosnoc:model:Base:p_global_val_size',...
                            'nosnoc: User provided p_global_val has the wrong size.')
                    end
                else
                    p_global_val = zeros(dims.n_p_global,1);
                end
            else
                dims.n_p_global = 0;
                obj.p_global = define_casadi_symbolic(opts.casadi_symbolic_mode,'',0);
                obj.p_global_val = [];
            end

            if size(obj.p_time_var, 1) ~= 0
                dims.n_p_time_var = size(obj.p_time_var, 1);
                if size(obj.p_time_var_val, 1) ~= 0
                    if size(obj.p_time_var_val) ~= [dims.n_p_time_var, opts.N_stages]
                        error('nosnoc:model:Base:p_time_var_val_size',...
                            'nosnoc: User provided p_time_var_val has the wrong size.')
                    end
                else
                    obj.p_time_var_val = zeros(dims.n_p_time_var, opts.N_stages);
                end
            else
                dims.n_p_time_var = 0;
                obj.p_time_var = define_casadi_symbolic(opts.casadi_symbolic_mode,'',0);
                obj.p_time_var_val = double.empty(0,opts.N_stages);
            end

            obj.p = vertcat(obj.p_global,obj.p_time_var);

            %% Stage and terminal costs check
            if size(obj.f_q, 1) == 0
                obj.f_q = 0;
            end

            if size(obj.f_q_T, 1) == 0
                obj.f_q_T = 0;
            end
            %% Least squares objective terms with variables references
            if size(obj.lsq_x, 1) ~= 0
                if length(obj.lsq_x)<3
                    error('nosnoc:model:Base:lsq_x_missing_component',...
                        'nosnoc: In lsq_x either the least squares function, the reference of the weight matrix are missing.')
                end
                if size(obj.lsq_x{2},1)~=size(obj.lsq_x{1})
                    error('nosnoc:model:Base:lsq_x_dim_mismatch',...
                        'nosnoc: The dimensions of the least squares error term and weighting matrix for the differential states do not match.')
                end
                if size(obj.lsq_x{1},1)~=size(obj.lsq_x{3})
                    error('nosnoc:model:Base:lsq_x_dim_mismatch',...
                        'nosnoc: The dimensions of the least squares error term and reference for the differential states do not match.')
                end

                n_x_ref_rows = size(obj.lsq_x{2},1);
                n_x_ref_cols = size(obj.lsq_x{2},2);
                if n_x_ref_cols == opts.N_stages
                    fprintf('nosnoc: the provided reference for the differential states is time variable. \n');
                elseif n_x_ref_cols == 1
                    % replaciate
                    fprintf('nosnoc: the provided reference for the differential states is constant over time. \n');
                    obj.lsq_x{2} = repmat(obj.lsq_x{2},1,opts.N_stages);
                else
                    fprintf('nosnoc: The reference in lsq_x has to have a length of %d (if constant) or %d if time variables. \n',1,opts.N_stages)
                    error('nosnoc:model:Base:lsq_x_dim_mismatch',...
                        'nosnoc: Please provide x_ref in lsq_x{1} with an appropriate size.')
                end
                obj.x_ref_val = obj.lsq_x{2};
                obj.x_ref = define_casadi_symbolic(opts.casadi_symbolic_mode,'x_ref',n_x_ref_rows);
                obj.f_lsq_x = (obj.lsq_x{1}-obj.x_ref)'*obj.lsq_x{3}*(obj.lsq_x{1}-obj.x_ref);
            else
                obj.x_ref = [];
                obj.f_lsq_x = [];
                obj.x_ref_val = [];
            end

            % least square terms for control inputs
            if size(obj.lsq_u, 1) ~= 0
                if length(obj.lsq_u)<3
                    error('nosnoc:model:Base:lsq_u_missing_component',...
                        'nosnoc: In lsq_u either the least squares function, the reference of the weight matrix are missing.')
                end
                if size(obj.lsq_u{2},1)~=size(obj.lsq_u{1})
                    error('nosnoc:model:Base:lsq_u_dim_mismatch',...
                        'nosnoc: The dimensions of the least squares error term and weighting matrix for the control input do not match.')
                end
                if size(obj.lsq_u{1},1)~=size(obj.lsq_u{3})
                    error('nosnoc:model:Base:lsq_u_dim_mismatch',...
                        'nosnoc: The dimensions of the least squares error term and reference for the control input do not match.')
                end
                n_u_ref_rows = size(obj.lsq_u{2},1);
                n_u_ref_cols = size(obj.lsq_u{2},2);
                if n_u_ref_cols == opts.N_stages
                    fprintf('nosnoc: the provided reference for the control inputs is time variable. \n');
                elseif n_u_ref_cols == 1
                    % replaciate
                    fprintf('nosnoc: the provided reference for the control inputs is constant over time. \n');
                    obj.lsq_u{2} = repmat(obj.lsq_u{2},1,opts.N_stages);
                else
                    fprintf('nosnoc: The reference in lsq_u has to have a length of %d (if constant) or %d if time variables. \n',1,opts.N_stages)
                    error('nosnoc:model:Base:lsq_u_dim_mismatch',...
                        'nosnoc: Please provide u_ref in lsq_u{2} with an appropriate size.')
                end
                obj.u_ref_val = obj.lsq_u{2};
                obj.u_ref = define_casadi_symbolic(opts.casadi_symbolic_mode,'u_ref',n_u_ref_rows);
                obj.f_lsq_u = (obj.lsq_u{1}-obj.u_ref)'*obj.lsq_u{3}*(obj.lsq_u{1}-obj.u_ref);
            else
                obj.u_ref = [];
                obj.f_lsq_u = [];
                obj.u_ref_val = [];
            end

            % least square terms for control inputs
            if size(obj.lsq_T, 1) ~= 0
                % sanity chkecs on the input
                if length(obj.lsq_T)<3
                    error('nosnoc:model:Base:lsq_T_missing_component',...
                        'nosnoc: In lsq_T either the least squares function, the reference or the weight matrix are missing.')
                end
                if size(obj.lsq_T{2},1)~=size(obj.lsq_T{1})
                    error('nosnoc:model:Base:lsq_T_dim_mismatch',...
                        'nosnoc: The dimensions of the least squares error term and weighting matrix for the terminal cost do not match.')
                end
                if size(obj.lsq_T{1},1)~=size(obj.lsq_T{3})
                    error('nosnoc:model:Base:lsq_T_dim_mismatch',...
                        'nosnoc: The dimensions of the least squares error term and reference for the terminal cost do not match.')
                end

                n_x_T_rows = size(obj.lsq_T{2},1);
                n_x_T_cols = size(obj.lsq_T{2},2);
                if n_x_T_cols == 1
                    fprintf('nosnoc: The provided reference for the terminal cost is ok. \n');
                else
                    fprintf('nosnoc: The reference in lsq_T has to be a vector of length %d. \n',length(obj.lsq_T{1}));
                    error('nosnoc:model:Base:lsq_T_dim_mismatch',...
                        'nosnoc: Please provide a reference vector in lsq_T{2} with an appropriate size.')
                end
                obj.x_ref_end_val = obj.lsq_T{2};
                obj.x_ref_end = define_casadi_symbolic(opts.casadi_symbolic_mode,'x_ref_end',n_x_T_rows);
                obj.f_lsq_T = (obj.lsq_T{1}-obj.x_ref_end)'*obj.lsq_T{3}*(obj.lsq_T{1}-obj.x_ref_end);
            else
                obj.x_ref_end = [];
                obj.f_lsq_T = [];
                obj.x_ref_end_val = [];
            end

            %% Inequality constraints check
            if size(obj.g_path, 1) ~= 0
                g_path_constraint  = 1;
                dims.n_g_path = length(obj.g_path);
                if size(obj.lbg_path, 1) ~= 0
                    if length(obj.lbg_path)~=dims.n_g_path;
                        error('nosnoc:model:Base:lbg_path_size',...
                            'The user provided vector lbg_path has the wrong size.')
                    end
                else
                    obj.lbg_path = -inf*ones(dims.n_g_path,1);
                end

                if size(obj.ubg_path, 1) ~= 0
                    if length(obj.ubg_path)~=dims.n_g_path;
                        error('nosnoc:model:Base:ubg_path_size',...
                            'The user provided vector ubg_path has the wrong size.')
                    end
                else
                    obj.ubg_path =  0*ones(dims.n_g_path,1);
                end
            else
                dims.n_g_path = 0;
                g_path_constraint = 0;
            end

            %% Check path complementarity constraints
            if size(obj.G_path, 1) ~= 0
                if size(obj.G_path, 1) ~= size(obj.H_path, 1)
                    error('nosnoc:model:Base:path_complementarity_mismatch',...
                        'G_path and H_path must be the same size.')
                end
            end
            %% Terminal constraints
            if size(obj.g_terminal, 1) ~= 0
                dims.n_g_terminal = length(obj.g_terminal);
                if size(obj.lbg_terminal, 1) ~= 0
                    if length(obj.lbg_terminal)~=dims.n_g_terminal
                        error('nosnoc:model:Base:lbg_terminal_size',...
                            'nosnoc: The provided vector lbg_terminal has the wrong size.')
                    end
                else
                    obj.lbg_terminal = 0*ones(dims.n_g_terminal,1);
                end

                if size(obj.ubg_terminal, 1) ~= 0
                    if length(obj.ubg_terminal)~=dims.n_g_terminal
                        error('nosnoc:model:Base:ubg_terminal_size',...
                            'nosnoc: The provided vector ubg_terminal has the wrong size.')
                    end
                else
                    obj.ubg_terminal =  0*ones(dims.n_g_terminal,1);
                end
            else
                dims.n_g_terminal = 0;
            end

            obj.dims = dims;
        end
    end
    % Print output in a structured way.
    methods(Access=protected)
        function propgrp = getPropertyGroups(obj)
            gTitle1 = 'Populated Properties';
            propList1 = struct;
            names = properties(obj);
            for ii=1:length(names)
                name = names{ii};
                
                if any(size(obj.(name)) == 0)
                    continue
                end

                % some custom handling for objective functions:
                if strcmp(name, 'f_q')
                    if class(obj.f_q) == "double" && obj.f_q == 0
                        continue
                    end
                end
                if strcmp(name, 'f_q_T')
                    if class(obj.f_q_T) == "double" && obj.f_q_T == 0
                        continue
                    end
                end
                
                propList1.(names{ii}) = obj.(name);% TODO(@anton) better custom display here
            end
            propgrp(1) = matlab.mixin.util.PropertyGroup(propList1,gTitle1);
        end

        function displayScalarObject(obj)
            className = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
            scalarHeader = [className ' Model'];
            header = sprintf('%s\n',scalarHeader);
            disp(header)
            propgroup = getPropertyGroups(obj);
            matlab.mixin.CustomDisplay.displayPropertyGroups(obj,propgroup)
        end
    end
end
