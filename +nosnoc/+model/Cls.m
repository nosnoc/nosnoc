classdef Cls < nosnoc.model.Base
% A class of systems of rigid bodies with contacts and friction with model equations:
%
% .. math::
%     :nowrap:
%
%     \begin{align*}
%        \dot{q}       &= v, \\
%        M\dot{v}     &= f_v(q,v) + \sum_{i=1}^{n_c} J_n^i \lambda_n^i + J_t^i \lambda_t^i, \\
%        0            &\le \lambda_n^i \perp f_c^i(q) \ge 0\\
%        0            &= J_n^i(q(t_s))^\top (v(t_s^+) + ev(t_s^-))\quad \mathrm{if}\ f_c^i(t_s) = 0\ \mathrm{and}\ J_n^i(q(t_s))^\top v(t_s^-) < 0 \\
%        \lambda_t^i  &\in \operatorname{arg\, min}_{\lambda_t^i\in\mathbb{R}^{n_t}} -v^\top  J_t^i \lambda_t^i\quad \mathrm{s.t.}\quad ||\lambda_t^i||_2 \le \mu^i \lambda_n^i
%     \end{align*}
%
% With $i = 1\ldots n_c$.
%
% See Also:
%     More information about this model and its discretization can be found in :cite:`Nurkanovic2024`.
    properties
        q % casadi.SX|casadi.MX: Gemeralized coordinates $q\in\mathbb{R}^{n_q}$.
        v % casadi.SX|casadi.MX: Gemeralized velocities $v\in\mathbb{R}^{n_q}$.
        f_v % casadi.SX|casadi.MX: Generalized acceleration $\dot{v} = f_v(x)\in\mathbb{R}^{n_q}$.
        f_c % casadi.SX|casadi.MX: Contact gap functions $f_c(q)\in\mathbb{R}^{n_c}$.

        % double: Friction coefficients $0\le\mu$. If a scalar value is passed then it is assumed all
        % contacts have the same friction coefficient. Otherwise this needs to be a vector of size $n_c$.
        mu

        % double: Coefficient of restituition $0\le e \le 1$. If a scalar value is passed then it is
        % assumed that all contact have the same coefficient of restitution. Otherwise this needs to
        % be a vector of size $n_c$.
        e

        M % casadi.SX|casadi.MX|double: Generalized inertia matrix. Can be a function of the state $q$.
        invM % casadi.SX|casadi.MX|double: Inverse of the generalized inertial matrix. TODO(@anton) should this be private/readonly

        friction_exists % boolean: Set to true if any $\mu \neq 0$.

        J_normal % casadi.SX|casadi.MX: Normal contact Jacobian $J_n$. This can be calculated automatically from the contact gap functions.
        J_tangent % casadi.SX|casadi.MX: Tangent contact Jacobian $J_t$. This must be provided if there is friction and using the Conic :attr:`~nosnoc.Options.friction_model`.

        % casadi.SX|casadi.MX: Tangent contact polyherdal approximation.
        % This must be provided if there is friction and using the Polyhedral :mat:attr:`~nosnoc.Options.friction_model`.
        % For every row $D_i$, $-D_i$ must also be a row in $D$.
        D_tangent 
    end

    methods
        function obj = Cls()
            obj = obj@nosnoc.model.Base();
        end
        
        function verify_and_backfill(obj, opts)
            arguments
                obj
                opts nosnoc.Options
            end
            import casadi.*
            verify_and_backfill@nosnoc.model.Base(obj,opts);

            dims = obj.dims;

            % TODO(@anton) generalizing this is a massive pain. Choose a strict interface.
            if isempty(obj.q)
                obj.q = obj.x(1:(dims.n_x/2));
            end

            if isempty(obj.v)
                obj.v = obj.x((dims.n_x/2 + 1):end);
            end

            dims.n_q = dims.n_x/2;
            dims.n_v = dims.n_x/2;
            
            if size(obj.f_v,1) ~= dims.n_v
                error("nosnoc: f_v has incorrect dimension. It must have the same dimension as x.")
            end

            dims.n_c = size(obj.f_c,1);

            if isempty(obj.e)
                error("nosnoc: Please provide a coefficient of restitution via model.e")
            else
                if length(obj.e) ~= 1 && length(obj.e) ~= dims.n_c
                    error('The length of model.e has to be one or match the length of model.f_c')
                end
                if length(obj.e) == 1
                    obj.e = obj.e*ones(dims.n_c,1);
                end
            end
            if any(abs(1-obj.e)>1) || any(obj.e<0)
                error('nosnoc: the coefficient of restitution e should be in [0,1].')
            end

            % coefficient of friction checks
            if size(obj.mu, 1) ~= 0
                if length(obj.mu) ~= 1 && length(obj.mu) ~= dims.n_c
                    error('The length of model.mu has to be one or match the length of model.f_c')
                end
                if length(obj.mu) == 1
                    obj.mu = obj.mu*ones(dims.n_c,1);
                end

                if any(obj.mu > 0)
                    obj.friction_exists = 1;
                else
                    obj.friction_exists = 0;
                end
            else
                obj.friction_exists = 0;
                obj.mu = zeros(dims.n_c,1);
                fprintf('nosnoc: Coefficients of friction mu not provided, setting it to zero for all contacts. \n')
            end
            if any(obj.mu<0)
                error('nosnoc: The coefficients of friction mu should be nonnegative.')
            end

            if isempty(obj.M)
                fprintf('nosnoc: Inertia matrix M is not provided, set to default value: M = eye(n_q). \n')
                obj.M = eye(dims.n_q);
                obj.invM = inv(obj.M);
            elseif any(size(obj.M) ~= dims.n_q)
                error('nosnoc: Inertia matrix has incorrect dimensions.')
            else
                obj.invM = inv(obj.M);
            end

            if size(obj.J_normal, 1) ~= 0
                J_normal_exists = 1;
            else
                J_normal_exists = 0;
            end

            if J_normal_exists
                if size(obj.J_normal,1)~=dims.n_q && size(obj.J_normal,2)~=dims.n_c
                    fprintf('nosnoc: J_normal should be %d x %d matrix.\n',dims.n_q,dims.n_c);
                    error('nosnoc: J_normal has the wrong size.')
                end
                J_normal_exists = 1;
            else
                obj.J_normal = obj.f_c.jacobian(obj.q)';
                fprintf('nosnoc: normal contact Jacobian not provided, but it is computed from the gap functions.\n');
                J_normal_exists = 1;
            end

            % Tangent Contact Jacobian
            if obj.friction_exists
                if isequal(opts.friction_model,'Conic')
                    if size(obj.J_tangent, 1) ~= 0
                        if size(obj.J_tangent,1)~=dims.n_q
                            error('nosnoc: J_tangent has the wrong size.')
                        end
                    else
                        error('nosnoc: please provide the tangent Jacobian in model.J_tangent.')
                    end
                end

                if isequal(opts.friction_model,'Polyhedral')
                    if isempty(obj.D_tangent)
                        error('nosnoc: please provide the polyhedral tangent Jacobian in model.D_tangent, e.g., using the conic tangent Jacobian model.J_tangent: D_tangent = [J_tangent(q_0),-J_tangent(q_0)].')
                    end
                end
            end
            % Dimension of tangents
            dims.n_t = 0;
            if obj.friction_exists
                if isequal(opts.friction_model,'Polyhedral')
                    dims.n_t = size(obj.D_tangent,2)/dims.n_c; % number of tangent multipliers for a single contact
                elseif isequal(opts.friction_model,'Conic')
                    dims.n_t = size(obj.J_tangent,2)/dims.n_c; % number of tangent multipliers for a single contact
                end
                dims.n_tangents = dims.n_t*dims.n_c; % number tangent forces for all multpliers
            else
                dims.n_tangents = 0;
            end

            obj.dims = dims;
        end
    end
end
