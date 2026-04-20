classdef Qpec < handle & matlab.mixin.Scalar
    % & matlab.mixin.CustomDisplay
    % Base class for Quadratic Programs with Equilibirumn Constraints (QPEC)
    %
    % .. math::
    %     :nowrap:
    %
    % \begin{equation*}
    %     \begin{aligned}
    %         \min_{w} \quad & w^\top Q^k w + w^\top w \\
    %         \text{s.t.} \quad &  g_{\mathrm{lb}} \leq b + \nabla g(x^k)^\top d \leq  g_{\mathrm{lb}},\\
    %         & 0\leq G(x^k) + \nabla G(x^k)^\top d \perp  G(x^k) + \nabla H(x^k)^\top d \geq 0.
    %     \end{aligned}
    % \end{equation*}

    % Which when comming from the linearization of MPCC looks like::
    % .. math::
    %     :nowrap:
    %
    % \begin{equation*}
    %     \begin{aligned}
    %         \min_{d} \quad & \frac{1}{2} d^\top H^k d + \nabla f(x^k)^\top d \\
    %         \text{s.t.} \quad &  g_{\mathrm{lb}} \leq g(x^k) + \nabla g(x^k)^\top d \leq  g_{\mathrm{lb}},\\
    %         & 0\leq G(x^k) + \nabla G(x^k)^\top d \perp  G(x^k) + \nabla H(x^k)^\top d \geq 0.
    %     \end{aligned}
    % \end{equation*}


    %
    % Represents and solves a linearized MPCC subproblem of the form:
    %
    %   minimize    0.5 * d' * Q * d + q' * d
    %   subject to  lba <= b + A*d <= uba,
    %               0 <= g + G*d ⟂ h + H*d >= 0,
    %               lbw <= d <= ubw.
    %
    % The class stores matrix data, bounds, and solver backend
    % (reg_homotopy, mpecopt, lcqpow, or gurobi).

    properties
        % ---- Sparsity patterns ----
        Q_sparsity   % Sparsity pattern of the quadratic term Q.
        A_sparsity   % Sparsity pattern of the linear constraint matrix A.
        G_sparsity   % Sparsity pattern of the complementarity constraint matrix G.
        H_sparsity   % Sparsity pattern of the complementarity constraint matrix H.

        % ---- Matrices ----
        Q            % Matrix for the quadratic cost term in the QPEC objective.
        A            % Jacobian of regular (inequality/equality) constraints.
        G            % Jacobian of the first complementarity function.
        H            % Jacobian of the second complementarity function.

        % ---- Vectors ----
        q            % Linear term of the quadratic objective.
        lbw          % Lower bound on decision variables w.
        ubw          % Upper bound on decision variables w.
        b            % Offset vector in regular constraints.
        lba          % Lower bounds on linear constraints A*w + b.
        uba          % Upper bounds on linear constraints A*w + b.
        g            % Offset vector in first complementarity constraint.
        h            % Offset vector in second complementarity constraint.

        % ---- Initialization ----
        w0           % Initial guess for QPEC decision variable.
        y0           % Integer variables if QPEC is solved as miqp or via sos1.

        % ---- Meta data and solver ----
        dims         % Struct containing problem dimensions (n_x, n_a, n_comp).
        qpec_solver       % Handle or object of the backend QPEC solver.
        qp_solver    % Handle or object of the backend QP solver.

        solver_opts  % Options structure for the selected QPEC solver.
        plugin       % Name of the solver backend (e.g., 'reg_homotopy', 'mpecopt', 'lcqpow', 'gurobi').

        qp_solver_opts  % Options structure for the selected QP solver.
        qp_plugin       % Name of the solver backend (e.g., qp_gurobi', 'qp_ipopt', 'qp_highs').

        % ---- Results and statistics ----
        result_qpec  % Result structure of the last QPEC solve.
        stats        % Statistics structure returned from the last solver call.
    end

    methods
        function obj = Qpec(Q, A, G, H)
            % Constructor: get sparsity patterns for QPEC matrices (needed mostly for Lcqpow)
            % NOTE: If passing numeric values for the matrices and they have a sparsity pattern: pass them as a sparse matrix,
            %       as otherwise we assume that they are dense which may lead to memory issues.
            import casadi.*
            % Process hessian sparsity
            if isa(Q, 'casadi.Sparsity')
                obj.Q_sparsity = Q;
                obj.Q = sparse(DM(obj.Q_sparsity));
            else
                obj.Q_sparsity = DM(Q).sparsity;
                obj.Q = Q;
            end
            % Process linear constraint sparsity
            if isa(A, 'casadi.Sparsity')
                obj.A_sparsity = A;
                obj.A = sparse(DM(obj.A_sparsity));
            else
                obj.A_sparsity = DM(A).sparsity;
                obj.A = A;
            end
            % Process complementarity constraint sparsity.
            if isa(G, 'casadi.Sparsity')
                obj.G_sparsity = G;
                obj.G = sparse(DM(obj.G_sparsity));
            else
                obj.G_sparsity = DM(G).sparsity;
                obj.G = G;
            end
            if isa(H, 'casadi.Sparsity')
                obj.H_sparsity = H;
                obj.H = sparse(DM(obj.H_sparsity));
            else
                obj.H_sparsity = DM(H).sparsity;
                obj.H = H;
            end

            % Get dimensions
            [dims.n_a, dims.n_x] = size(obj.A_sparsity);
            dims.n_comp = size(obj.G_sparsity, 1);
            obj.dims = dims;

            % TODO(@anton) check sizes

            % Populate vector values (allocate dummy values)
            obj.q = zeros(dims.n_x,1);
            obj.lbw = -inf(dims.n_x,1);
            obj.ubw = inf(dims.n_x,1);
            obj.w0 = zeros(dims.n_x,1);
            obj.lba = zeros(dims.n_a,1);
            obj.uba = inf(dims.n_a,1);
            obj.g = zeros(dims.n_comp,1);
            obj.h = zeros(dims.n_comp,1);
            obj.y0 = zeros(dims.n_comp,1);
        end

        function create_qpec_solver(obj, opts, plugin)
            % Create and configure a QPEC solver for the selected backend.
            arguments
                obj
                opts
                plugin = 'reg_homotopy'
            end
            obj.plugin = plugin;
            obj.solver_opts = opts;

            switch plugin
                case {'reg_homotopy', 'mpecopt', 'ccopt'}
                    mpec = obj.create_parametric_qpec;
                    obj.qpec_solver = nosnoc.solver.mpccsol('qpec_solver', plugin, mpec, opts);
                case {'lcqpow'}
                    % Lcqpow solver is re-initialized every call.

                case {'gurobi'}
                    % Build Gurobi QPCC solver backend
                    obj.qpec_solver = nosnoc.qpec.GurobiQpecSolver(opts);
                otherwise
                    nosnoc.error('invalid_qpec_plugin', ...
                        ['Plugin ' plugin ' is not valid. Use: reg_homotopy, mpecopt, lcqpow, or gurobi.'])
            end
        end

        function create_qp_solver(obj, opts, plugin)
            % Create and configure a QPEC solver for the selected backend.
            arguments
                obj
                opts
                plugin = 'qp_ipopt'
            end
            obj.qp_plugin = plugin;
            obj.qp_solver_opts = opts;

            switch plugin
                case {'qp_gurobi'}
                    % Build Gurobi QPCC solver backend
                    obj.qp_solver = nosnoc.qpec.GurobiQpecSolver(opts);
                    % todo: port Rtmpc qpsolver creation here;
                case {'qp_ipopt'}
                    % assumes that ipopt options are passed to qpec solver constructor
                    import casadi.*
                    qp = obj.create_parametric_qp();
                    obj.qp_solver = nlpsol('solver', 'ipopt', qp, opts);
                case {'qp_highs'}
                    error('not implemented yet.')
                otherwise
                    obj.qp_solver = [];
            end
        end

        function fname = codegen_qpec(obj, fname)
            import casadi.*
            % Generate data
            dims = obj.dims;
            n_a = dims.n_a;
            n_comp = dims.n_comp;
            n_a_aug = n_a + 2*n_comp;
            x = SX.sym('x', dims.n_x);
            p = SX.sym('p', 0);
            lam_g = SX.sym('lam_g', n_a_aug);
            lam_f = SX.sym('lam_f', 1);
            f = 0.5*x'*obj.Q*x + obj.q'*x;
            g = obj.A*x + obj.b;
            G = obj.G*x + obj.g;
            H = obj.H*x + obj.h;
            g_aug = vertcat(g, G, H);
            grad_f = f.gradient(x);
            jac_g = g_aug.jacobian(x);
            L = lam_f*f - lam_g'*g_aug;
            [hess_L, nabla_L] = L.hessian(x);

            % Generate funs
            nlp_f = Function('nlp_f', {x, p}, {f}, {'x','p'}, {'f'});
            nlp_grad_f = Function('nlp_grad_f', {x, p}, {f, grad_f}, {'x','p'}, {'f','grad_f'});
            nlp_g = Function('nlp_g', {x, p}, {g_aug}, {'x','p'}, {'g'});
            nlp_jac_g = Function('nlp_jac_g', {x, p}, {g, jac_g}, {'x', 'p'}, {'g', 'jac_g'});
            nlp_hess_l = Function('nlp_hess_l', {x, p, lam_f, lam_g}, {hess_L}, {'x', 'g', 'lam_f', 'lam_g'}, {'hess_l'});

            % Generate json
            json_struct.x0 = obj.w0;
            json_struct.y0 = zeros(n_a_aug, 1);
            json_struct.p0 = [];
            json_struct.lbx = obj.lbw;
            json_struct.ubx = obj.ubw;
            json_struct.lbg = vertcat(obj.lba, zeros(2*n_comp, 1));
            json_struct.ubg = vertcat(obj.uba, inf*ones(2*n_comp, 1));
            json_struct.ind_cc1 = n_a+1:n_a+n_comp;
            json_struct.ind_cc2 = n_a+n_comp+1:n_a+2*n_comp;
            json = jsonencode(json_struct, "ConvertInfAndNaN", false, "PrettyPrint", true);
            fid = fopen([fname, '.json'], "w");
            fprintf(fid, json);
            fclose(fid);

            % Generate c code
            cg = CodeGenerator([fname,'.c']);
            cg.add(nlp_f);
            cg.add(nlp_g);
            cg.add(nlp_grad_f);
            cg.add(nlp_jac_g);
            cg.add(nlp_hess_l);
            cg.generate();
        end

        function [result_qpec, stats] = solve(obj)
            dims = obj.dims;
            % Solve the current QPEC using the chosen solver backend.
            switch obj.plugin
                case {'reg_homotopy', 'mpecopt', 'ccopt'}
                    qpec_initialization.lbg = obj.lba;
                    qpec_initialization.ubg = obj.uba;
                    qpec_initialization.lbx = obj.lbw;
                    qpec_initialization.ubx = obj.ubw;
                    qpec_initialization.x0 = obj.w0; % TODO maybe have initialization a object param
                    % NOTE: all QPEC data is stored in a parameter, this avoids creating a new solver when using mpccsol/mpecopt which have ipopt underneath
                    qpec_initialization.p = vertcat(obj.Q(:), obj.q, obj.A(:), obj.G(:), obj.H(:), obj.b, obj.g, obj.h);
                    result_qpec = obj.qpec_solver(qpec_initialization);
                    stats = obj.qpec_solver.stats;
                    %  stats.homotopy_iterations todo : use this as nice QPCC warm start
                case {'lcqpow'}
                    lba = obj.lba - obj.b;
                    uba = obj.uba - obj.b;
                    [w_res, lam_res, stats] = LCQPow(sparse(obj.Q), obj.q, sparse(obj.G), sparse(obj.H), -obj.g, inf(dims.n_comp,1), -obj.h, inf(dims.n_comp,1), sparse(obj.A), lba, uba, obj.lbw, obj.ubw, obj.solver_opts);
                    result_qpec.x = w_res;
                    % TODO(@anton) this is a guess!
                    result_qpec.f = 0.5*w_res'*obj.Q*w_res + obj.q'*w_res;
                    result_qpec.lam_x = lam_res(1:dims.n_x);
                    result_qpec.lam_g = lam_res(dims.n_x + (1:dims.n_a));
                case {'gurobi'}
                    % Solve via Gurobi backend (MIQP/SOS1/REG)
                    [result_qpec, stats] = obj.qpec_solver.solve(obj);
                otherwise
                    nosnoc.error('invalid_qpec_plugin', ['Plugin ' plugin ' is not valid please use: reg_homotopy, mpecopt, lcqpow, or gurobi'])
            end
            obj.result_qpec = result_qpec;
            obj.stats = stats;
        end

        function [result_qp, stats] = solve_qp(obj)
            dims = obj.dims;
            % Solve the current QPEC using the chosen solver backend.
            switch obj.qp_plugin
                case {'qp_ipopt'}
                    p_val = vertcat(obj.Q(:), obj.q, obj.A(:), obj.G(:), obj.H(:), obj.b, obj.g, obj.h);
                    % p_val = sparse(vertcat(obj.Q(:), obj.q, obj.A(:), obj.G(:), obj.H(:), obj.b, obj.g, obj.h, obj.y0));
                    % result_qp = obj.qp_solver(qp_initialization);
                    lbg_qp = [obj.lba; zeros(2*dims.n_comp,1)];
                    ubg_qp = obj.uba;
                    for ii = 1:dims.n_comp
                        if obj.y0(ii) == 1
                            ubg_qp  = [ubg_qp; inf];
                            ubg_qp  = [ubg_qp; 0];
                        else
                            ubg_qp  = [ubg_qp; 0];
                            ubg_qp  = [ubg_qp; inf];
                        end
                    end
                    % result_qp = obj.qp_solver('lbx',obj.lbw,'ubx',obj.ubw,'lbg',lbg_qp,'ubg',ubg_qp,'x0',obj.w0,'p',p_val);
                    result_qp = obj.qp_solver('lbx',obj.lbw,'ubx',obj.ubw,'lbg',lbg_qp,'ubg',ubg_qp,'p',p_val);
                    stats = obj.qp_solver.stats();
                case {'qp_highs'}
                    error('not implemented yet.')
                case {'qp_gurobi'}
                    % Solve QP via Gurobi
                    [result_qp, stats] = obj.qp_solver.solve(obj);
                otherwise
                    nosnoc.error('invalid_qpec_plugin', ['Plugin ' plugin ' is not valid please use: qp_ipopt, qp_highs, or qp_gurobi'])
            end
            obj.result_qpec = result_qp;
            obj.stats = stats;
        end
    end

    methods(Access=private)
        function qpec = create_parametric_qpec(obj)
            % Needed only for mpccsol and mpecopt!

            % This functions creates a mpcc (a qpec/qpcc, where all data are casadi parameters).
            % This qoec struct is passed to mpccsol or mpecopt, and the values of the parameters are updated every call.
            % This way we create a mpccsol/mpecopt solver only once.

            % Form of QPEC:
            % min 0.5w'* Q_p* w + q_p' w
            % s.t.  lba <= A_p*w  + b_p <= uba
            %       lbw <= w <= ubw;
            %        0 <= g_p + G_p*w \perp h_p + h_p*w >=0

            import casadi.*
            dims = obj.dims;

            q_p = SX.sym('q_p', dims.n_x, 1);
            Q_p = SX.sym('Q_p', obj.Q_sparsity);
            % Regular constraint data parameters
            b_p = SX.sym('b_p', dims.n_a, 1);
            A_p = SX.sym('A_p', obj.A_sparsity);
            % Complementarity constraints data parameters
            g_p = SX.sym('G_p', dims.n_comp, 1);
            G_p = SX.sym('nabla_G_p', obj.G_sparsity);
            h_p = SX.sym('H_p', dims.n_comp, 1);
            H_p = SX.sym('nabla_H_p', obj.H_sparsity);

            p_qpec = vertcat(Q_p(:),q_p,A_p(:),G_p(:),H_p(:),b_p,g_p,h_p);

            w = SX.sym('w', dims.n_x); % qpec decision variables

            % store data
            qpec = struct;
            qpec.x = w;
            qpec.f = 0.5*w'*Q_p*w + q_p'*w;%+rho_d*w'*w;;
            qpec.g = b_p + A_p*w;
            qpec.G = g_p + G_p*w;
            qpec.H = h_p + H_p*w;
            qpec.p = p_qpec; % all problem parameters
        end

        function qp = create_parametric_qp(obj)
            % Needed for parametric branch QP if y is given

            % This functions creates a mpcc (a qpec/qpcc, where all data are casadi parameters).
            % This qoec struct is passed to mpccsol or mpecopt, and the values of the parameters are updated every call.
            % This way we create a mpccsol/mpecopt solver only once.

            % Form of QPEC:
            % min 0.5w'* Q_p* w + q_p' w
            % s.t.  lba <= A_p*w  + b_p <= uba
            %       lbw <= w <= ubw;
            %        0 <= g_p + G_p*w \perp h_p + h_p*w >=0

            import casadi.*
            dims = obj.dims;

            q_p = SX.sym('q_p', dims.n_x, 1);
            Q_p = SX.sym('Q_p', obj.Q_sparsity);
            % Regular constraint data parameters
            b_p = SX.sym('b_p', dims.n_a, 1);
            A_p = SX.sym('A_p', obj.A_sparsity);
            % Complementarity constraints data parameters
            g_p = SX.sym('G_p', dims.n_comp, 1);
            G_p = SX.sym('nabla_G_p', obj.G_sparsity);
            h_p = SX.sym('H_p', dims.n_comp, 1);
            H_p = SX.sym('nabla_H_p', obj.H_sparsity);
            p_qpec = vertcat(Q_p(:),q_p,A_p(:),G_p(:),H_p(:),b_p,g_p,h_p);

            w = SX.sym('w', dims.n_x); % qpec decision variables

            % Represent both as comps with using y:
            g_expr = g_p + G_p*w;
            h_expr = h_p + H_p*w;
            GH_comp = [];
            for ii = 1:dims.n_comp
                GH_comp = [GH_comp; g_expr(ii);h_expr(ii)];
            end

            % store data
            qp = struct;
            qp.x = w;
            qp.f = 0.5*w'*Q_p*w + q_p'*w;%+rho_d*w'*w;;
            qp.g = [b_p + A_p*w; GH_comp];
            qp.p = p_qpec; % all problem parameters
        end
    end
end
