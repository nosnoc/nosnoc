% BSD 2-Clause License

% Copyright (c) 2023, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% This file is part of NOSNOC.
classdef Objects < Base
properties (Access=public)
        % Differential state
        x
        lbx
        ubx
        x0

        % user algebraics
        z
        z0
        lbz
        ubz
        g_z % user algebraic constraints

        % Controls
        u
        lbu
        ubu
        u0

        % Algebraics
        lambda
        mu
        alpha
        p_d
        y1_d
        y2_d
        lambda_t
        v_t
        normal_lift
        gamma_f

        % x_dot lift
        x_dot_lift

        % friction expressions
        G_friction
        H_friction
        g_friction

        % dynamics
        f_x

        % Feasible set
        r = 0.01
        g_d
        g_kkt
        normal
        normal0
        tangent
        tangent0

        % global variables
        v_global
        v0_global
        lbv_global
        ubv_global
        
        % global parameters
        p_global
        p_global_val

        % time varying parameters
        p_time_var
        p_time_var_val
        % params
        p_time_var_stages
        % All params
        p

        % TODO: maybe a separate OCP class common to all model types?
        % Objective
        f_q
        f_q_T

        % least squares
        lsq_x
        x_ref
        f_lsq_x
        x_ref_val
        lsq_u
        u_ref
        f_lsq_u
        u_ref_val
        lsq_T
        x_ref_end
        f_lsq_T
        x_ref_end_val

        % Path constraints
        g_path
        lbg_path
        ubg_path

        % Terminal constraints
        g_terminal
        lbg_terminal
        ubg_terminal

        % Dimensions
        dims

        % Distance functions
        c

        % Functions
        f_x_fun
        c_fun
        g_kkt_fun
        g_d_fun
        g_tangent_fun
        g_friction_fun
        G_friction_fun
        H_friction_fun
        normal_fun
        g_z_fun
        g_path_fun
        g_terminal_fun
        f_q_fun
        f_q_T_fun

        % objects
        objects
        contacts

        %
        verified = 0
        variables_created = 0
        functions_created = 0
    end

    methods (Access=public)
        function obj = Objects()
            dims = struct;
        end

        function generate_functions(obj)
            import casadi.*
            if obj.functions_created
                return
            end
            dims = obj.dims;

            obj.f_x = [];

            for ii=1:length(obj.objects)
                obj.x_dot_lift = vertcat(obj.x_dot_lift, obj.objects(ii).x_dot_lift);
                obj.f_x = vertcat(obj.f_x, obj.objects(ii).x_dot);
            end

            if 1 % TODO(@anton) add option for lifting x_dot
                opts.allow_free = true;
                g_friction_sub_fun = Function('g_friction_sub', {obj.x_dot_lift}, {obj.g_friction}, opts);
                obj.g_friction = g_friction_sub_fun(obj.f_x);
            end

            obj.g_kkt_fun = Function('g_kkt', {obj.x, obj.mu, obj.p_d, obj.y1_d, obj.y2_d}, {obj.g_kkt});
            obj.g_d_fun = Function('g_d', {obj.x, obj.alpha, obj.p_d, obj.y1_d, obj.y2_d}, {obj.g_d});
            obj.normal_fun = Function('normal', {obj.x, obj.mu, obj.p_d, obj.p}, {obj.normal});
            obj.g_friction_fun = Function('g_friction', {obj.x, obj.u, obj.lambda, obj.lambda_t, obj.mu, obj.p_d, obj.normal_lift, obj.tangent, obj.v_t, obj.p}, {obj.g_friction});
            obj.G_friction_fun = Function('G_friction', {obj.lambda, obj.lambda_t, obj.v_t, obj.gamma_f}, {obj.G_friction});
            obj.H_friction_fun = Function('H_friction', {obj.lambda_t, obj.gamma_f}, {obj.H_friction});
            obj.f_x_fun = Function('f_x', {obj.x, obj.z, obj.u, obj.lambda, obj.lambda_t, obj.tangent, obj.mu, obj.p_d, obj.v_global, obj.p}, {obj.f_x});
            obj.c_fun = Function('c', {obj.x, obj.alpha, obj.v_global, obj.p}, {obj.c});
            obj.g_z_fun = Function('g_z', {obj.x, obj.z, obj.u, obj.v_global, obj.p}, {obj.g_z});
            obj.g_path_fun = Function('g_path', {obj.x, obj.z, obj.u, obj.v_global, obj.p_global, obj.p_time_var}, {obj.g_path});
            obj.g_terminal_fun = Function('g_terminal', {obj.x, obj.z, obj.v_global, obj.p_global}, {obj.g_terminal});
            obj.f_q_fun = Function('f_q', {obj.x, obj.z, obj.u, obj.v_global, obj.p}, {obj.f_q});
            obj.f_q_T_fun = Function('f_q_T', {obj.x, obj.z, obj.v_global, obj.p}, {obj.f_q_T});

            obj.functions_created = true;
        end

        function generate_variables(obj)
            if obj.variables_created
                return
            end
            obj.variables_created = true;
        end

        function verify_and_backfill(obj, opts)
            import casadi.*
            if obj.verified
                return
            end

            dims = obj.dims;

            %% Check differential state
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

                % check x0
                if size(obj.x0, 1) ~= 0
                    if length(obj.x0) ~= dims.n_x
                        error('nosnoc: The vector x0, for the initial value of x has the wrong size.')
                    end
                else
                    obj.x0 = zeros(dims.n_x,1);
                end
            else
                error('nosnoc: Please provide the state vector x, a CasADi symbolic variable.');
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
                obj.z = define_casadi_symbolic(opts.casadi_symbolic_mode,'',0);
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
                obj.u = define_casadi_symbolic(opts.casadi_symbolic_mode,'',0);
                obj.u0 = [];
                dims.n_u = 0;
                if opts.print_level >=1
                    fprintf('nosnoc: No control vector u iif isa(shape1, 'nosnoc.objects.Polygon') && isa(shape2, 'nosnoc.objects.Polygon')
                obj.addPolygonPolygon(shape1, shape2, mu);
            elseif isa(shape1, 'nosnoc.objects.Polygon') && isa(shape2, 'nosnoc.objects.Ball')
                obj.addPolygonBall(shape1, shape2, mu);
            elseif isa(shape1, 'nosnoc.objects.Ball') && isa(shape2, 'nosnoc.objects.Polygon')
                obj.addPolygonBall(shape2, shape1, mu);
            elses provided. \n')
                end
                obj.lbu = [];
                obj.ubu = [];
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
                obj.p_global = define_casadi_symbolic(opts.casadi_symbolic_mode,'',0);
                obj.p_global_val = [];
            end

            if size(obj.p_time_var, 1) ~= 0
                dims.n_p_time_var = size(obj.p_time_var, 1);
                if size(obj.p_time_var_val, 1) ~= 0
                    if size(obj.p_time_var_val) ~= [dims.n_p_time_var, opts.N_stages]
                        error('nosnoc: User provided p_time_var_val has the wrong size.')
                    end
                else
                    obj.p_time_var_val = zeros(dims.n_p_time_var, opts.N_stages);
                end

                obj.p_time_var_stages = [];
                for ii=1:opts.N_stages
                    var_full = define_casadi_symbolic(opts.casadi_symbolic_mode, ['p_time_var_' num2str(ii)], dims.n_p_time_var);
                    obj.p_time_var_stages = horzcat(obj.p_time_var_stages, var_full);
                end
            else
                dims.n_p_time_var = 0;
                obj.p_time_var = define_casadi_symbolic(opts.casadi_symbolic_mode,'',0);
                obj.p_time_var_stages = define_casadi_symbolic(opts.casadi_symbolic_mode,'', [0, opts.N_stages]);
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
                    error('nosnoc: In lsq_x either the least squares function, the reference of the weight matrix are missing.')
                end
                if size(obj.lsq_x{2},1)~=size(obj.lsq_x{1})
                    error('nosnoc: The dimensions of the least squares error term and weighting matrix for the differential states do not match.')
                end
                if size(obj.lsq_x{1},1)~=size(obj.lsq_x{3})
                    error('nosnoc: The dimensions of the least squares error term and reference for the differential states do not match.')
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
                    fprintf('nosnoc: The reference in lsq_x has to have a length of %d (if constant) or %d if time vriables. \n',1,opts.N_stages)
                    error('nosnoc: Please provide x_ref in lsq_x{1} with an appropriate size.')
                end
                obj.x_ref_val = obj.lsq_x{2};
                obj.x_ref = define_casadi_symbolic(opts.casadi_symbolic_mode,'x_ref',n_x_ref_rows);
                obj.f_lsq_x = (obj.lsq_x{1}-obj.x_ref)'*obj.lsq_x{3}*(obj.lsq_x{1}-obj.x_ref);
            else
                obj.x_ref = define_casadi_symbolic(opts.casadi_symbolic_mode,'x_ref',1);
                obj.f_lsq_x = 0;
                obj.x_ref_val = zeros(1,opts.N_stages);
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
                if n_u_ref_cols == opts.N_stages
                    fprintf('nosnoc: the provided reference for the control inputs is time variable. \n');
                elseif n_u_ref_cols == 1
                    % replaciate
                    fprintf('nosnoc: the provided reference for the control inputs is constant over time. \n');
                    obj.lsq_u{2} = repmat(obj.lsq_u{2},1,opts.N_stages);
                else
                    fprintf('nosnoc: The reference in lsq_u has to have a length of %d (if constant) or %d if time vriables. \n',1,opts.N_stages)
                    error('nosnoc: Please provide u_ref in lsq_u{2} with an appropriate size.')
                end
                obj.u_ref_val = obj.lsq_u{2};
                obj.u_ref = define_casadi_symbolic(opts.casadi_symbolic_mode,'u_ref',n_u_ref_rows);
                obj.f_lsq_u = (obj.lsq_u{1}-obj.u_ref)'*obj.lsq_u{3}*(obj.lsq_u{1}-obj.u_ref);
            else
                obj.u_ref = define_casadi_symbolic(opts.casadi_symbolic_mode,'u_ref',1);
                obj.f_lsq_u = 0;
                obj.u_ref_val = zeros(1,opts.N_stages);
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

                n_x_T_rows = size(obj.lsq_T{2},1);
                n_x_T_cols = size(obj.lsq_T{2},2);
                if n_x_T_cols == 1
                    fprintf('nosnoc: the provided reference for the terminal cost is ok. \n');
                else
                    fprintf('nosnoc: The reference in lsq_T has to be a vector of length %d. \n',length(obj.lsq_T{1}));
                    error('nosnoc: Please provide a reference vector in lsq_T{2} with an appropriate size.')
                end
                obj.x_ref_end_val = obj.lsq_T{2};
                obj.x_ref_end = define_casadi_symbolic(opts.casadi_symbolic_mode,'x_ref_end',n_x_T_rows);
                obj.f_lsq_T = (obj.lsq_T{1}-obj.x_ref_end)'*obj.lsq_T{3}*(obj.lsq_T{1}-obj.x_ref_end);
            else
                obj.x_ref_end  = define_casadi_symbolic(opts.casadi_symbolic_mode,'x_ref_end',1);
                obj.f_lsq_T = 0;
                obj.x_ref_end_val = 0;
            end

            %% v_global
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
                    if length(obj.ubv_global) ~= n_v_global
                        error('nosnoc: The vector ubv_global, for the upper bound of v_global has the wrong size.')
                    end
                else
                    obj.ubv_global = inf*ones(n_v_global, 1);
                end
            else
                n_v_global = 0;
                obj.v_global = define_casadi_symbolic(opts.casadi_symbolic_mode, '', 0);
                obj.v0_global = [];
                obj.lbv_global = [];
                obj.ubv_global = [];
            end
            dims.n_v_global = n_v_global;

            %% Inequality constraints check
            if size(obj.g_path, 1) ~= 0
                g_path_constraint  = 1;
                n_g_path = length(obj.g_path);
                if size(obj.lbg_path, 1) ~= 0
                    if length(obj.lbg_path)~=n_g_path;
                        error('The user provided vector lbg_path has the wrong size.')
                    end
                else
                    obj.lbg_path = -inf*ones(n_g_path,1);
                end

                if size(obj.ubg_path, 1) ~= 0
                    if length(obj.ubg_path)~=n_g_path;
                        error('The user provided vector ubg_path has the wrong size.')
                    end
                else
                    obj.ubg_path =  0*ones(n_g_path,1);
                end
            else
                n_g_path = 0;
                g_path_constraint  = 0;
            end

            %% Terminal constraints
            if size(obj.g_terminal, 1) ~= 0
                n_g_terminal = length(obj.g_terminal);
                if size(obj.lbg_terminal, 1) ~= 0
                    if length(obj.lbg_terminal)~=n_g_terminal
                        error('nosnoc: The provided vector lbg_terminal has the wrong size.')
                    end
                else
                    obj.lbg_terminal = 0*ones(n_g_terminal,1);
                end

                if size(obj.ubg_terminal, 1) ~= 0
                    if length(obj.ubg_terminal)~=n_g_terminal
                        error('nosnoc: The provided vector ubg_terminal has the wrong size.')
                    end
                else
                    obj.ubg_terminal =  0*ones(n_g_terminal,1);
                end
            else
                n_g_terminal = 0;
            end

            dims.n_c = size(obj.c, 1);
            dims.n_y1d = size(obj.y1_d, 1);
            dims.n_y2d = size(obj.y2_d, 1);
            dims.n_mu = size(obj.mu, 1);
            dims.n_pd = size(obj.p_d,1);
            dims.n_alpha = size(obj.alpha,1);
            dims.n_tangent = size(obj.tangent,1);
            dims.n_normal = size(obj.normal,1);
            dims.n_normal_lift = size(obj.normal_lift,1);
            dims.n_v_t = size(obj.v_t,1);
            dims.n_gamma_f = size(obj.gamma_f,1);
            dims.n_lambda_t = size(obj.lambda_t,1);
            obj.dims = dims;

            obj.verified = true;
        end

        function addContact(obj, shape1, shape2, mu)
            arguments
                obj
                shape1 nosnoc.objects.Object
                shape2 nosnoc.objects.Object
                mu = 0
            end
            if shape1.n_dim ~= shape2.n_dim
                error('Dimensions mismatch');
            end
            
            if isa(shape1, 'nosnoc.objects.Ball') && isa(shape2, 'nosnoc.objects.Ball')
                obj.addBallBall(shape1, shape2, mu);
            elseif isa(shape1, 'nosnoc.objects.Ellipse') && isa(shape2, 'nosnoc.objects.Ellipse')
                obj.addEllipseEllipse(shape1, shape2, mu);
            elseif isa(shape1, 'nosnoc.objects.Ball') && isa(shape2, 'nosnoc.objects.Ellipse')
                obj.addEllipseBall(shape2,shape1, mu);
            elseif isa(shape1, 'nosnoc.objects.Ellipse') && isa(shape2, 'nosnoc.objects.Ball')
                obj.addEllipseBall(shape1,shape2, mu);
            else
                error('Unknown shape types');
            end

            if ~ismember(shape1, obj.objects)
                obj.objects = [obj.objects, shape1];
            end

            if ~ismember(shape2, obj.objects)
                obj.objects = [obj.objects, shape2];
            end

            obj.contacts{end+1} = {shape1, shape2, mu};

            obj.x = [];
            obj.x0 = [];
            obj.lbx = [];
            obj.ubx = [];
            for ii=1:length(obj.objects)
                obj.x = vertcat(obj.x, obj.objects(ii).x);
                obj.x0 = vertcat(obj.x0, obj.objects(ii).x0);
                obj.lbx = vertcat(obj.lbx, obj.objects(ii).lbx);
                obj.ubx = vertcat(obj.ubx, obj.objects(ii).ubx);
            end
        end
    end

    methods (Access=private)
        function addBallBall(obj, ball1, ball2, mu_f)
            import casadi.*

            % Lambda and gap functions
            lambda = SX.sym(['lambda_' ball1.name '_' ball2.name]);
            c = sum((ball1.x - ball2.x).^2) - (ball1.r + ball2.r)^2;

            % Get normals
            nabla_c_x1 = c.jacobian(ball1.x)';
            nabla_c_x2 = c.jacobian(ball2.x)';

            % Update dynamics of both balls
            ball1.x_dot = ball1.x_dot + nabla_c_x1*lambda;
            ball2.x_dot = ball2.x_dot + nabla_c_x2*lambda;

            normal = [nabla_c_x1;nabla_c_x2];
            obj.normal = [obj.normal;normal];
            obj.normal0 = [obj.normal0;0;1;0;-1];

            if mu_f
                if ball1.n_dim == 2
                    tangent = SX.sym('tangent', 2);
                    obj.tangent = [obj.tangent; tangent];
                    obj.tangent0 = [obj.tangent0;1;0];
                    v_t = SX.sym('v_t', 2);
                    obj.v_t = [obj.v_t; v_t];
                    lambda_t = SX.sym('lambda_t', 2);
                    obj.lambda_t = [obj.lambda_t; lambda_t];
                    gamma_f = SX.sym('gamma_f', 1);
                    obj.gamma_f = [obj.gamma_f; gamma_f];
                    normal_lift = SX.sym('normal_lift', 4);
                    obj.normal_lift = [obj.normal_lift; normal_lift];

                    ball1.x_dot = ball1.x_dot + lambda_t(1)*tangent - lambda_t(2)*tangent;
                    ball2.x_dot = ball2.x_dot - lambda_t(1)*tangent + lambda_t(2)*tangent;
                    
                    g_friction = [normal_lift - normal; % Lift normal (helps with convergence) TODO(@anton) allow for not lifting
                        tangent - [-normal_lift(2);normal_lift(1)]; % get tangent (lifted again) TODO(@anton) allow for not lifting
                        v_t(1) - (dot(ball1.x_dot_lift, tangent)/dot(tangent,tangent)); % Tangent velocities of both objects 
                        v_t(2) - (dot(ball2.x_dot_lift, tangent)/dot(tangent,tangent))];
                    obj.g_friction = [obj.g_friction; g_friction];
                    
                    G_friction = [v_t(1) - v_t(2) + gamma_f;
                        v_t(2) - v_t(1) + gamma_f;
                        mu_f*lambda - lambda_t(1) - lambda_t(2);
                        lambda_t(1)]; % TODO(@anton) can we avoid this as it violates MPCC-LICQ
                    obj.G_friction = [obj.G_friction; G_friction];
                    H_friction = [lambda_t; gamma_f;lambda_t(2)];
                    obj.H_friction = [obj.H_friction; H_friction];
                    
                else
                    error("3d friction not yet supported")
                end
            end
            
            % update model
            obj.c = vertcat(obj.c,c);
            obj.lambda = vertcat(obj.lambda, lambda);
        end

        function addEllipseBall(obj, ellipse, ball, mu_f)
            import casadi.*
            n_dim = ball.n_dim;

            % Helper Functions
            if n_dim == 2
                theta = SX.sym('theta');
                R_matrix = [cos(theta) -sin(theta);...
                    sin(theta) cos(theta)];
                R = Function('R', {theta}, {R_matrix});
                a = SX.sym('a',2);
                b = SX.sym('b',2);
                cross_fun = Function('cross', {a,b}, {a(1)*b(2) - a(2)*b(1)});
            else
                rx = SX.sym('rx',1);
                ry = SX.sym('ry',1);
                rz = SX.sym('rz',1);

                Rot = [cos(ry)*cos(rz), sin(rx)*sin(ry)*cos(rz) - cos(rx)*sin(rz),  cos(rx)*sin(ry)*cos(rz) + sin(rx)*sin(rz);
                    cos(ry)*sin(rz), sin(rx)*sin(ry)*sin(rz) + cos(rx)*cos(rz),  cos(rx)*sin(ry)*sin(rz) - sin(rx)*cos(rz);
                    -sin(ry), sin(rx)*cos(ry),  cos(rx)*cos(ry)];
                R = Function('R', {[rx;ry;rz]}, {Rot});
                a = SX.sym('a',3);
                b = SX.sym('b',3);
                cross_fun = Function('cross', {a,b}, {cross(a,b)});
            end

            % Define contacts
            for ii=1:ellipse.N
                lambda = SX.sym(['lambda_' ellipse.name '_' ball.name]);
                obj.lambda = vertcat(obj.lambda, lambda);
                alpha = SX.sym(['alpha_' ellipse.name '_' ball.name]);
                obj.alpha = vertcat(obj.alpha, alpha);
                p_d = SX.sym(['p_d_' ellipse.name '_' ball.name], n_dim);
                obj.p_d = vertcat(obj.p_d, p_d);
                mu = SX.sym(['mu_' ellipse.name '_' ball.name], 2);
                obj.mu = vertcat(obj.mu, mu);

                A = ellipse.A{ii};
                r = ball.r;

                mu_1c = mu(end-1);
                mu_2c = mu(end);
                
                g_d = [(p_d-ellipse.c)'*R(ellipse.xi)*A*R(ellipse.xi)'*(p_d-ellipse.c) - alpha;
                      (1/r^2)*(p_d-ball.c)'*(p_d-ball.c) - alpha];
                obj.g_d = vertcat(obj.g_d, g_d);

                L = alpha - 1 + mu'*g_d;
                
                % g_kkt = [1 - mu_1c - mu_2c;
                %     mu_1c*(R(ellipse.xi)*(A + A')*R(ellipse.xi)')*(p_d-ellipse.c) + mu_2c*(2/r^2)*(p_d-ball.c)];
                g_kkt = L.jacobian([alpha;p_d])';
                obj.g_kkt = vertcat(obj.g_kkt, g_kkt);

                
                ntr = 2/ball.r^2*(p_d - ball.c)*mu_2c;
                normal_ball = [-ntr];
                normal_ellipse =[ntr;
                    cross_fun(p_d-ellipse.c, ntr)];

                obj.normal = [obj.normal;normal_ball;normal_ellipse];
                obj.normal0 = [obj.normal0;0;1;0;-1;0];

                % Update dynamics of both objects
                ellipse.x_dot = ellipse.x_dot + normal_ellipse*lambda;
                ball.x_dot = ball.x_dot + normal_ball*lambda;

                if mu_f
                    if ball.n_dim == 2 %TODO does this need changes?
                        tangent = SX.sym('tangent', 2);
                        obj.tangent = [obj.tangent; tangent];
                        obj.tangent0 = [obj.tangent0;1;0];
                        v_t = SX.sym('v_t', 2);
                        obj.v_t = [obj.v_t; v_t];
                        lambda_t = SX.sym('lambda_t', 2);
                        obj.lambda_t = [obj.lambda_t; lambda_t];
                        gamma_f = SX.sym('gamma_f', 1);
                        obj.gamma_f = [obj.gamma_f; gamma_f];
                        normal_lift = SX.sym('normal_lift', 4);
                        obj.normal_lift = [obj.normal_lift; normal_lift];
                        
                        ellipse.x_dot = ellipse.x_dot + lambda_t(1)*[tangent;cross_fun(p_d-ellipse.c, tangent)] - lambda_t(2)*[tangent;cross_fun(p_d-ellipse.c, tangent)];
                        ball.x_dot = ball.x_dot - lambda_t(1)*tangent + lambda_t(2)*tangent;
                        
                        g_friction = [normal_lift - [normal_ball;normal_ellipse(1:2)]; % Lift normal (helps with convergence) TODO(@anton) allow for not lifting
                            tangent - [-normal_lift(2);normal_lift(1)]; % get tangent (lifted again) TODO(@anton) allow for not lifting
                            v_t(1) - (dot(ellipse.x_dot_lift(1:2), tangent)/dot(tangent,tangent) + ellipse.x_dot_lift(3)*norm(p_d - ellipse.c)); % Tangent velocities of both objects 
                            v_t(2) - (dot(ball.x_dot_lift, tangent)/dot(tangent,tangent))];
                        obj.g_friction = [obj.g_friction; g_friction];
                        
                        G_friction = [v_t(1) - v_t(2) + gamma_f;
                            v_t(2) - v_t(1) + gamma_f;
                            mu_f*lambda - lambda_t(1) - lambda_t(2)]; % TODO(@anton) can we avoid this as it violates MPCC-LICQ
                        obj.G_friction = [obj.G_friction; G_friction];
                        H_friction = [lambda_t; gamma_f];
                        obj.H_friction = [obj.H_friction; H_friction];
                    else
                        error("3d friction not yet supported")
                    end
                end
                
                obj.c = vertcat(obj.c, alpha - 1);
            end
        end
        
        function addEllipseEllipse(obj, ellipse1, ellipse2, mu_f)
            import casadi.*
            n_dim = ellipse2.n_dim;

            % Helper Functions
            if n_dim == 2
                theta = SX.sym('theta');
                R_matrix = [cos(theta) -sin(theta);...
                    sin(theta) cos(theta)];
                R = Function('R', {theta}, {R_matrix});
                a = SX.sym('a',2);
                b = SX.sym('b',2);
                cross_fun = Function('cross', {a,b}, {a(1)*b(2) - a(2)*b(1)});
            else
                rx = SX.sym('rx',1);
                ry = SX.sym('ry',1);
                rz = SX.sym('rz',1);

                Rot = [cos(ry)*cos(rz), sin(rx)*sin(ry)*cos(rz) - cos(rx)*sin(rz),  cos(rx)*sin(ry)*cos(rz) + sin(rx)*sin(rz);
                    cos(ry)*sin(rz), sin(rx)*sin(ry)*sin(rz) + cos(rx)*cos(rz),  cos(rx)*sin(ry)*sin(rz) - sin(rx)*cos(rz);
                    -sin(ry), sin(rx)*cos(ry),  cos(rx)*cos(ry)];
                R = Function('R', {[rx;ry;rz]}, {Rot});
                a = SX.sym('a',3);
                b = SX.sym('b',3);
                cross_fun = Function('cross', {a,b}, {cross(a,b)});
            end

            % Define contacts
            for ii=1:ellipse1.N
                for jj=1:ellipse2.N
                    lambda = SX.sym(['lambda_' ellipse1.name '_' ellipse2.name]);
                    obj.lambda = vertcat(obj.lambda, lambda);
                    alpha = SX.sym(['alpha_' ellipse1.name '_' ellipse2.name]);
                    obj.alpha = vertcat(obj.alpha, alpha);
                    p_d = SX.sym(['p_d_' ellipse1.name '_' ellipse2.name], n_dim);
                    obj.p_d = vertcat(obj.p_d, p_d);
                    mu = SX.sym(['mu_' ellipse1.name '_' ellipse2.name], 2);
                    obj.mu = vertcat(obj.mu, mu);

                    A1 = ellipse1.A{ii};
                    A2 = ellipse2.A{jj};

                    mu_1c = mu(end-1);
                    mu_2c = mu(end);
                    
                    g_d = [(p_d-ellipse1.c)'*R(ellipse1.xi)*A1*R(ellipse1.xi)'*(p_d-ellipse1.c) - alpha;
                        (p_d-ellipse2.c)'*R(ellipse2.xi)*A2*R(ellipse2.xi)'*(p_d-ellipse2.c) - alpha];
                    obj.g_d = vertcat(obj.g_d, g_d);

                    L = alpha - 1 + mu'*g_d;
                    
                    % g_kkt = [1 - mu_1c - mu_2c;
                    %     mu_1c*(R(ellipse1.xi)*(A + A')*R(ellipse1.xi)')*(p_d-ellipse1.c) + mu_2c*(2/r^2)*(p_d-ellipse2.c)];
                    g_kkt = L.jacobian([alpha;p_d])';
                    obj.g_kkt = vertcat(obj.g_kkt, g_kkt);

                    
                    ntr = -2*R(ellipse1.xi)*A1*R(ellipse1.xi)'*(p_d - ellipse1.c)*mu_1c;
                    normal_ellipse1 =[ntr;
                        cross_fun(p_d-ellipse1.c, ntr)];
                    normal_ellipse2 = [-ntr;
                        -cross_fun(p_d-ellipse2.c, ntr)];

                    obj.normal = [obj.normal;normal_ellipse1;normal_ellipse2];
                    obj.normal0 = [obj.normal0;0;1;0;-1;0];

                    % Update dynamics of both objects
                    ellipse1.x_dot = ellipse1.x_dot + normal_ellipse1*lambda;
                    ellipse2.x_dot = ellipse2.x_dot + normal_ellipse2*lambda;

                    if mu_f
                        if ellipse2.n_dim == 2 %TODO does this need changes?
                            tangent = SX.sym('tangent', 2);
                            obj.tangent = [obj.tangent; tangent];
                            obj.tangent0 = [obj.tangent0;1;0];
                            v_t = SX.sym('v_t', 2);
                            obj.v_t = [obj.v_t; v_t];
                            lambda_t = SX.sym('lambda_t', 2);
                            obj.lambda_t = [obj.lambda_t; lambda_t];
                            gamma_f = SX.sym('gamma_f', 1);
                            obj.gamma_f = [obj.gamma_f; gamma_f];
                            normal_lift = SX.sym('normal_lift', 4);
                            obj.normal_lift = [obj.normal_lift; normal_lift];
                            
                            ellipse1.x_dot = ellipse1.x_dot + lambda_t(1)*[tangent;cross_fun(p_d-ellipse1.c, tangent)] - lambda_t(2)*[tangent;cross_fun(p_d-ellipse1.c, tangent)];
                            ellipse2.x_dot = ellipse2.x_dot - lambda_t(1)*tangent + lambda_t(2)*tangent;
                            
                            g_friction = [normal_lift - [normal_ellipse2;normal_ellipse1(1:2)]; % Lift normal (helps with convergence) TODO(@anton) allow for not lifting
                                tangent - [-normal_lift(2);normal_lift(1)]; % get tangent (lifted again) TODO(@anton) allow for not lifting
                                v_t(1) - (dot(ellipse1.x_dot_lift(1:2), tangent)/dot(tangent,tangent) + ellipse1.x_dot_lift(3)*norm(p_d - ellipse1.c)); % Tangent velocities of both objects 
                                v_t(2) - (dot(ellipse2.x_dot_lift, tangent)/dot(tangent,tangent))];
                            obj.g_friction = [obj.g_friction; g_friction];
                            
                            G_friction = [v_t(1) - v_t(2) + gamma_f;
                                v_t(2) - v_t(1) + gamma_f;
                                mu_f*lambda - lambda_t(1) - lambda_t(2);
                                lambda_t(1);
                                lambda]; % TODO(@anton) can we avoid this as it violates MPCC-LICQ
                            obj.G_friction = [obj.G_friction; G_friction];
                            H_friction = [lambda_t; gamma_f;lambda_t(2);gamma_f];
                            obj.H_friction = [obj.H_friction; H_friction];
                        else
                            error("3d friction not yet supported")
                        end
                    end
                    
                    obj.c = vertcat(obj.c, alpha - 1);
                end
            end
        end
    end
end
