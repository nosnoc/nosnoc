classdef GurobiQpecSolver < handle
    % Gurobi-based solver for QPCCs (Quadratic Programs with Complementarity Constraints).
    % Stateless: each call to `solve(qpec)` builds and solves a new Gurobi model.
    %
    % Supports methods: 'miqp' (Big-M), 'sos1', 'reg' (Scholtes regularization), and 'bilinear'

    properties
        options    % qpec.GurobiOptions instance
        result     % struct with x, y, lam_g, lam_x, objective, etc.
        stats      % solver info
    end

    methods
        function obj = GurobiQpecSolver(options)
            % Solver constructor
            if nargin < 1 || isempty(options)
                options = qpec.GurobiOptions(); %
            end
            obj.options = options;
        end

        function [result, stats] = solve(obj, qpec)
            % Top-level dispatcher
            t_solve = tic;
            switch obj.options.method
                case {'qp'}
                    % solves a branch QP based on the given y (active set);
                    [result, stats] = obj.solve_qp_fixed_y(qpec);
                case {'miqp'}
                    [result, stats] = obj.solve_miqp(qpec);
                case 'sos1'
                    [result, stats] = obj.solve_sos1(qpec);
                case 'reg'
                    [result, stats] = obj.solve_reg(qpec);
                case 'bilinear'
                    error('GurobiQpecSolver: Bilinear formulation not implemented yet.');
                otherwise
                    error('Unknown Gurobi method: %s', obj.options.method);
            end
            stats.cpu_time = toc(t_solve);
            obj.result = result;
            obj.stats  = stats;
        end
    end

    methods (Access = private)
        %% =================================================================
        % MIQP (Big-M) : build + solve + parse
        function [result, stats] = solve_miqp(obj, qpec)
            n  = qpec.dims.n_x;
            m  = qpec.dims.n_a;
            nc = qpec.dims.n_comp;

            % Big-M
            M = obj.options.bigM * ones(nc,1);

            % --- Model ---
            model = struct();
            model.modelsense = 'min';
            model.vtype = [repmat('C', n, 1); repmat('B', nc, 1)];
            model.lb    = [qpec.lbw; zeros(nc,1)];
            model.ub    = [qpec.ubw; ones(nc,1)];
            model.Q     = blkdiag(sparse(qpec.Q), sparse(nc,nc));
            model.obj   = [qpec.q; zeros(nc,1)];

            % --- Constraints (preserve row order) ---
            Arows = {}; rhs = {}; sense = {};
            A=qpec.A; b=qpec.b; lba=qpec.lba; uba=qpec.uba;

            % Regular rows
            for i = 1:m
                bi = b(i);
                if isfinite(lba(i)) && isfinite(uba(i)) && lba(i)==uba(i)
                    Arows{end+1,1} = [ A(i,:), sparse(1,nc) ];
                    rhs  {end+1,1} =  uba(i) - bi;
                    sense{end+1,1} = '=';
                else
                    if isfinite(uba(i))
                        Arows{end+1,1} = [ A(i,:), sparse(1,nc) ];
                        rhs  {end+1,1} =  uba(i) - bi;
                        sense{end+1,1} = '<';
                    end
                    if isfinite(lba(i))
                        Arows{end+1,1} = [-A(i,:), sparse(1,nc) ];
                        rhs  {end+1,1} = -(lba(i) - bi);
                        sense{end+1,1} = '<';
                    end
                end
            end

            % Big-M complementarity blocks
            if nc > 0
                G=qpec.G; H=qpec.H; g=qpec.g; h=qpec.h;
                % g + G w <= M y   ->  G w - M y <= -g
                Arows{end+1,1} = [ G, -spdiags(M,0,nc,nc) ];
                rhs  {end+1,1} = -g;
                sense{end+1,1} = repmat('<',nc,1);
                % h + H w <= M(1-y) -> H w + M y <= M - h
                Arows{end+1,1} = [ H,  spdiags(M,0,nc,nc) ];
                rhs  {end+1,1} = (M - h);
                sense{end+1,1} = repmat('<',nc,1);

                % nonnegativity for comp sides: g+Gw >=0 -> -G w <= g ; h+Hw >=0 -> -H w <= h
                if obj.options.lower_bounds_comps
                    Arows{end+1,1} = [-G, sparse(nc,nc)]; rhs{end+1,1} =  g; sense{end+1,1} = repmat('<', nc,1); %#ok<AGROW>
                    Arows{end+1,1} = [-H, sparse(nc,nc)]; rhs{end+1,1} =  h; sense{end+1,1} = repmat('<', nc,1); %#ok<AGROW>
                end
            end

            model.A     = sparse(vertcat(Arows{:}));
            model.rhs   = vertcat(rhs{:});
            model.sense = vertcat(sense{:});

            % --- Warm start (w0,y0) ---
            if obj.options.warmstart
                model.start(1:n) = qpec.w0;
                model.start(n+(1:nc)) = qpec.y0;
            end

            % --- Solve ---
            params = obj.options.gurobi_params;
            grb_result = gurobi(model, params);

            % --- Parse primal ---
            result = obj.parse_result(grb_result, qpec);

            % --- Stats (standardized shape) ---
            stats = struct();
            stats.status        = obj.field_or(grb_result,'status','unknown');
            stats.objval        = obj.field_or(grb_result,'objval',NaN);
            stats.objbound      = obj.field_or(grb_result,'objbound',NaN);
            stats.mipgap        = obj.field_or(grb_result,'mipgap',NaN);
            stats.runtime       = obj.field_or(grb_result,'runtime',NaN);
            stats.nodecount     = obj.field_or(grb_result,'nodecount',NaN);
            stats.itercount     = obj.field_or(grb_result,'itercount',NaN);
            stats.baritercount  = obj.field_or(grb_result,'baritercount',NaN);
            stats.solcount      = obj.field_or(grb_result,'solcount',NaN);
            stats.numvars       = numel(model.lb);
            stats.numbinaries   = nc;
            stats.numconstr     = size(model.A,1);
            stats.numsos        = 0;
            stats.mipgap_rel    = abs((stats.objval - stats.objbound)/max(1,abs(stats.objval)));
            stats.gurobi_raw    = grb_result;
            stats.success       = isfield(grb_result,'status') && strcmpi(grb_result.status,'OPTIMAL');


            % --- Dual recovery (optional) ---
            result.lam_g = [];
            result.lam_x = [];
            if obj.options.recover_duals && obj.options.resolve_qp_for_duals
                if stats.success
                    qpec.w0 =  grb_result.x(1:n);
                    qpec.y0  = grb_result.x(n+1:end);
                    [qpresults, qpstats] = obj.solve_qp_fixed_y(qpec, params);
                    result.lam_g = qpresults.lam_g;
                    result.lam_x = qpresults.lam_x;
                    stats.qp_recovery = qpstats;
                    stats.dual_source = 'resolved_qp';
                end
            end

            % --- Log ---
            persistent miqp_log_header_printed
            if isempty(miqp_log_header_printed)
                fprintf('%-10s %-12s %12s %12s %12s %10s\n', 'Solve','Status','ObjVal','||x||_inf','Runtime','Nodes');
                miqp_log_header_printed = true;
            end
            if isfield(grb_result,'x')
                fprintf('%-10s %-12s %12.4g %12.4g %12.4g %10.0f\n', ...
                    'QPCC-MIQP', stats.status, stats.objval, norm(result.x,inf), stats.runtime, stats.nodecount);
            else
                fprintf('%-10s %-12s %12s %12s %12.4g %10.0f\n', ...
                    'QPCC-MIQP', stats.status, 'NaN', 'NaN', stats.runtime, stats.nodecount);
            end
        end

        %% =================================================================
        % SOS1 : build + solve + parse (dual recovery via fixed-y QP)
        function [result, stats] = solve_sos1(obj, qpec)
            n  = qpec.dims.n_x;
            m  = qpec.dims.n_a;
            nc = qpec.dims.n_comp;

            % Vars: [w ; sG ; sH], sG,sH >= 0, SOS1([sG_i,sH_i]) per i
            % sG and sH are slack variables;
            if obj.options.lower_bounds_comps
                lb = [qpec.lbw(:); zeros(nc,1); zeros(nc,1)];
            else
                lb = [qpec.lbw(:); -inf(nc,1); -inf(nc,1)];
            end
            ub = [qpec.ubw(:);  inf(nc,1);  inf(nc,1)];
            vtype = repmat('C', n + 2*nc, 1);

            % --- Constraints (preserve MIQP row pattern) ---
            Arows = {}; rhs = {}; sense = {};
            A=qpec.A; b=qpec.b; lba=qpec.lba; uba=qpec.uba;

            % Regular rows
            for i = 1:m
                bi = b(i);
                if isfinite(lba(i)) && isfinite(uba(i)) && lba(i)==uba(i)
                    Arows{end+1,1} = [ A(i,:),  sparse(1,2*nc) ];
                    rhs  {end+1,1} =  uba(i) - bi;
                    sense{end+1,1} = '=';
                else
                    if isfinite(uba(i))
                        Arows{end+1,1} = [ A(i,:),  sparse(1,2*nc) ];
                        rhs  {end+1,1} =  uba(i) - bi;
                        sense{end+1,1} = '<';
                    end
                    if isfinite(lba(i))
                        Arows{end+1,1} = [-A(i,:),  sparse(1,2*nc) ];
                        rhs  {end+1,1} = -(lba(i) - bi);
                        sense{end+1,1} = '<';
                    end
                end
            end

            % Complementarity linearization: G*w + g - sG = 0, H*w + h - sH = 0
            for i = 1:nc
                e = sparse(1,i,1,1,nc);
                Arows{end+1,1} = [qpec.G(i,:), -e, sparse(1,nc)];
                rhs{end+1,1} = -qpec.g(i); sense{end+1,1} = '=';
            end
            for i = 1:nc
                e = sparse(1,i,1,1,nc);
                Arows{end+1,1} = [qpec.H(i,:), sparse(1,nc), -e];
                rhs{end+1,1} = -qpec.h(i); sense{end+1,1} = '=';
            end

            % --- Model ---
            model = struct();
            model.modelsense = 'min';

            model.A     = sparse(vertcat(Arows{:}));
            model.rhs   = vertcat(rhs{:});
            model.sense = vertcat(sense{:});
            model.Q     = blkdiag(sparse(qpec.Q), sparse(2*nc,2*nc));
            model.obj   = [qpec.q(:); zeros(2*nc,1)];
            model.lb    = lb;
            model.ub    = ub;
            model.vtype = vtype;

            % SOS1 sets
            sos = repmat(struct('type',1,'index',[],'weight',[]), nc, 1);
            for i = 1:nc
                sos(i).index  = [n+i, n+nc+i];   % sG_i, sH_i
                sos(i).weight = [1,     2     ];
            end
            model.sos = sos;

            % Warm start
            if obj.options.warmstart
                start = zeros(n + 2*nc, 1);
                if ~isempty(qpec.w0)
                    start(1:n) = qpec.w0(:);
                end
                if ~isempty(qpec.y0)
                    sG0 = (qpec.y0(:) < 0.5);
                    sH0 = (qpec.y0(:) > 0.5);
                    % sG0 = max(0, qpec.g + qpec.G*qpec.w0);
                    % sH0 = max(0, qpec.h + qpec.H*qpec.w0);
                    start(n+(1:nc))        = sG0;
                    start(n+nc+(1:nc))     = sH0;
                end
                if any(~isnan(start))
                    model.start = start;
                end
            end

            % Params
            params = obj.options.gurobi_params;
            params.Cuts = 0;
            params.Heuristics = 0;
            % Solve
            t0 = tic;
            gr = gurobi(model, params);
            cpu = toc(t0);

            % Extract primal
            w  = zeros(n,1);
            sG = zeros(nc,1);
            sH = zeros(nc,1);
            if isfield(gr,'x')
                x  = gr.x(:);
                w  = x(1:n);
                sG = x(n+(1:nc));
                sH = x(n+nc+(1:nc));
            end
            % Derive y from slacks
            tol = 1e-8;
            y_star = double(sG <= sH - tol);

            % Dual recovery: fix y and solve continuous QP for multiplier
            if obj.options.recover_duals && obj.options.resolve_qp_for_duals && isfield(grb_result,'status') && strcmpi(grb_result.status,'OPTIMAL');
                [qpresult, qpstats] = obj.solve_qp_fixed_y(qpec, params);
                lam_x = qpresult.lam_x;
                lam_g = qpresult.lam_g;
            else
                lam_x = zeros(n,1);
                lam_g = zeros(m,1);
            end

            % Result
            result = struct('x', w, 'y', y_star, 'lam_g', lam_g, 'lam_x', lam_x, 'obj', obj.field_or(gr,'objval',NaN));

            % Stats (same structure as MIQP)
            stats = struct();
            stats.status        = obj.field_or(gr,'status','unknown');
            stats.objval        = obj.field_or(gr,'objval',NaN);
            stats.objbound      = obj.field_or(gr,'objbound',NaN);   % may be empty for SOS, keep NaN fallback
            stats.mipgap        = obj.field_or(gr,'mipgap',NaN);     % may be empty
            stats.runtime       = obj.field_or(gr,'runtime',NaN);
            stats.nodecount     = obj.field_or(gr,'nodecount',NaN);
            stats.itercount     = obj.field_or(gr,'itercount',NaN);
            stats.baritercount  = obj.field_or(gr,'baritercount',NaN);
            stats.solcount      = obj.field_or(gr,'solcount',NaN);
            stats.numvars       = numel(model.lb);
            stats.numbinaries   = 0;               % no binaries in SOS1
            stats.numconstr     = size(model.A,1);
            stats.numsos        = numel(model.sos);
            stats.mipgap_rel    = abs((stats.objval - stats.objbound)/max(1,abs(stats.objval)));
            stats.gurobi_raw    = gr;
            stats.success       = strcmpi(stats.status,'OPTIMAL');
            if obj.options.recover_duals && obj.options.resolve_qp_for_duals
                stats.qp_recovery   = qpstats;
            end

            % Log
            persistent sos_log_header_printed
            if isempty(sos_log_header_printed)
                fprintf('%-10s %-12s %12s %12s %12s %10s\n', 'Solve','Status','ObjVal','||x||_inf','Runtime','Nodes');
                sos_log_header_printed = true;
            end
            if isfield(gr,'x')
                fprintf('%-10s %-12s %12.4g %12.4g %12.4g %10.0f\n', 'QPCC-SOS1', stats.status, stats.objval, norm(result.x,inf), stats.runtime, stats.nodecount);
            else
                fprintf('%-10s %-12s %12s %12s %12.4g %10.0f\n', 'QPCC-SOS1', stats.status, 'NaN','NaN', stats.runtime, stats.nodecount);
            end
        end

        %% =================================================================
        function [result, stats] = solve_reg(obj, qpec)
            % Scholtes regularization (SQP in absolute variable x):
            % minimize      0.5 x' Q x + q' x
            % s.t.          lba <= b + A x <= uba
            %               g+Gx >= 0,  h+Hx >= 0
            %               (g+Gx).*(h+Hx) <= tau   (linearized at current xk)
            %
            % Linearization at xk:
            %   let gG = g+G*xk,  gH = h+H*xk,  Jc = diag(gH)*G + diag(gG)*H
            %   then  (g+Gx).*(h+Hx) ≈ gG.*gH + Jc*(x - xk)  <= tau
            %   => Jc*x <= tau - gG.*gH + Jc*xk
            switch_after = 2; % when to swich to active set QP: TODO make option
            % ---- options ----
            tau        = obj.options.tau0;
            kappa      = obj.options.kappa;
            N_homotopy = obj.options.N_homotopy;
            Nsqp       = obj.options.Nsqp;
            sqp_tol    = obj.options.sqp_tol;
            verbose    = obj.options.verbose;

            % ---- dims & data ----
            n  = qpec.dims.n_x;
            m  = qpec.dims.n_a;
            nc = qpec.dims.n_comp;

            Q  = sparse(qpec.Q);              % symmetric (as you guarantee)
            q  = qpec.q(:);
            A  = qpec.A;  b = qpec.b;  lba = qpec.lba;  uba = qpec.uba;
            G  = qpec.G;  H = qpec.H;  g  = qpec.g;    h  = qpec.h;
            lb = qpec.lbw(:);  ub = qpec.ubw(:);

            % ---- Gurobi params (convex QPs) ----
            params = obj.options.gurobi_params;
            if ~isfield(params,'OutputFlag'),     params.OutputFlag = double(verbose); end
            if ~isfield(params,'Method'),         params.Method = 2;  end   % barrier
            if ~isfield(params,'Crossover'),      params.Crossover = 1; end % get pi/rc
            if ~isfield(params,'DualReductions'), params.DualReductions = 0; end
            if ~isfield(params,'NonConvex'),      params.NonConvex = 0; end

            params_bar = params;
            params_bar.Method = 2;  % barrier
            params_act = params;
            params_act.Method = 1;  % dual simplex


            if verbose
                fprintf('%-10s %4s %4s %12s %12s %12s %10s\n', ...
                    'Solve','ih','it','tau','||c||_inf','||Δx||_inf','Runtime');
            end

            % ---- initialize iterate xk ----
            xk = qpec.w0(:); if isempty(xk), xk = zeros(n,1); end

            t_all = tic;
            total_qp_calls = 0;
            total_qp_time  = 0;
            stages = repmat(struct(), N_homotopy, 1);

            last_gr = struct();
            last_map = struct('Aeq',nan(m,1),'Aub',nan(m,1),'Alb',nan(m,1));

            for ih = 1:N_homotopy
                t_stage = tic;
                last_obj = NaN; res_inf = Inf; step_inf = Inf;
                qptime_acc = 0; qp_count = 0;

                for it = 1:Nsqp
                    % ---- eval at xk ----
                    Axk = A*xk;
                    gG  = G*xk + g;
                    gH  = H*xk + h;

                    % product at xk (for monitoring)
                    c_xk = gG .* gH - tau;

                    % Jacobian of product at xk
                    Jc = spdiags(gH,0,nc,nc)*G + spdiags(gG,0,nc,nc)*H;

                    % ---- build QP in absolute x ----
                    model = struct();
                    model.modelsense = 'min';
                    model.Q   = Q;
                    model.obj = q;

                    Arows = {}; rhs = {}; sense = {};
                    row_map = struct('Aeq',nan(m,1),'Aub',nan(m,1),'Alb',nan(m,1));

                    % (1) ordinary rows: lba <= b + A x <= uba
                    for i = 1:m
                        ai = A(i,:); bi = b(i);
                        lb_i = lba(i); ub_i = uba(i);
                        has_lb = isfinite(lb_i); has_ub = isfinite(ub_i);

                        if has_lb && has_ub && lb_i == ub_i
                            Arows{end+1,1} =  ai;
                            rhs  {end+1,1} =  ub_i - bi;
                            sense{end+1,1} =  '=';
                            row_map.Aeq(i) = numel(rhs);
                        else
                            if has_ub
                                Arows{end+1,1} =  ai;
                                rhs  {end+1,1} =  ub_i - bi;
                                sense{end+1,1} =  '<';
                                row_map.Aub(i) = numel(rhs);
                            end
                            if has_lb
                                Arows{end+1,1} = -ai;
                                rhs  {end+1,1} = -(lb_i - bi);
                                sense{end+1,1} =  '<';
                                row_map.Alb(i) = numel(rhs);
                            end
                        end
                    end

                    % (2) nonnegativity of complementarity sides at x:
                    %     g+Gx >= 0  ->  -G x <= g
                    %     h+Hx >= 0  ->  -H x <= h
                    if nc > 0
                        Arows{end+1,1} = -G; rhs{end+1,1} =  g; sense{end+1,1} = repmat('<',nc,1);
                        Arows{end+1,1} = -H; rhs{end+1,1} =  h; sense{end+1,1} = repmat('<',nc,1);
                    end

                    % (3) Scholtes linearization at xk:
                    %     gG.*gH + Jc*(x - xk) <= tau  ->
                    %     Jc*x <= tau - gG.*gH + Jc*xk
                    if nc > 0
                        Arows{end+1,1} =  Jc;
                        rhs  {end+1,1} =  tau - (gG .* gH) + Jc*xk;
                        sense{end+1,1} =  repmat('<',nc,1);
                    end

                    % assemble
                    model.A     = sparse(vertcat(Arows{:}));
                    model.rhs   = vertcat(rhs{:});
                    model.sense = char(vertcat(sense{:}));

                    % bounds
                    model.lb = lb;
                    model.ub = ub;
                    model.start = xk;

                    % ---- solve QP ----
                    t_qp = tic;
                    if ih <= switch_after
                        gr = gurobi(model, params_bar);
                    else
                        model.start = xk;
                        gr = gurobi(model, params_act);
                    end

                    qp_t = toc(t_qp);
                    total_qp_calls = total_qp_calls + 1;
                    total_qp_time  = total_qp_time  + qp_t;
                    qptime_acc     = qptime_acc     + qp_t;

                    if ~isfield(gr,'x'), break; end
                    x = gr.x(:);

                    step_inf = norm(x - xk, inf);
                    gG_new   = G*x + g;
                    gH_new   = H*x + h;
                    res_inf  = max(0, max(gG_new .* gH_new - tau));
                    last_obj = obj.field_or(gr,'objval',NaN);
                    qp_count = qp_count + 1;

                    last_gr  = gr;
                    last_map = row_map;

                    if verbose
                        fprintf('%-10s %4d %4d %12.3g %12.3g %12.3g %10.3g\n', ...
                            'QPCC-Reg', ih, it, tau, res_inf, step_inf, qp_t);
                    end

                    xk = x;
                    if (res_inf <= max(tau, sqp_tol)) && (step_inf <= sqp_tol)
                        break;
                    end
                end

                stages(ih).tau         = tau;
                stages(ih).objval      = last_obj;
                stages(ih).runtime     = toc(t_stage);
                stages(ih).c_violation = res_inf;
                stages(ih).step_inf    = step_inf;
                stages(ih).qp_calls    = qp_count;

                tau = kappa * tau;
            end

            % ---- multipliers for original linear rows (optional) ----
            lam_g = zeros(m,1); lam_x = zeros(n,1);
            if isfield(last_gr,'pi')
                pi = last_gr.pi(:);
                for i = 1:m
                    if ~isnan(last_map.Aeq(i))
                        lam_g(i) = pi(last_map.Aeq(i));
                    else
                        ui = last_map.Aub(i); li = last_map.Alb(i);
                        if ~isnan(ui) && ~isnan(li)
                            lam_g(i) = pi(ui) - pi(li);
                        elseif ~isnan(ui)
                            lam_g(i) = pi(ui);
                        elseif ~isnan(li)
                            lam_g(i) = -pi(li);
                        end
                    end
                end
            end
            if isfield(last_gr,'rc'), lam_x = last_gr.rc(:); end

            % ---- heuristic y from final x ----
            y = [];
            if nc > 0
                gGf = G*xk + g;
                gHf = H*xk + h;
                tol = 1e-10;
                y = double(gHf < gGf - tol);
            end

            % ---- stats ----
            stats = struct();
            stats.status         = obj.field_or(last_gr,'status','unknown');
            stats.objval         = obj.field_or(last_gr,'objval',NaN);
            stats.runtime        = toc(t_all);
            stats.total_qp_calls = total_qp_calls;
            stats.total_qp_time  = total_qp_time;
            stats.per_stage      = stages;
            stats.method         = 'Scholtes (product<=tau, SQP in x)';
            stats.success        = isfield(last_gr,'status') && strcmpi(last_gr.status,'OPTIMAL');

            % ---- one-line log ----
            persistent reg_tail_hdr
            if isempty(reg_tail_hdr)
                fprintf('%-10s %-12s %12s %12s %12s %10s\n', ...
                    'Solve','Status','ObjVal','||x||_inf','Runtime','Nodes');
                reg_tail_hdr = true;
            end
            fprintf('%-10s %-12s %12.4g %12.4g %12.4g %10.0f\n', ...
                'QPCC-Reg', stats.status, stats.objval, norm(xk,inf), stats.runtime, total_qp_calls);

            % ---- result ----
            result = struct('x', xk, 'y', y, 'lam_g', lam_g, 'lam_x', lam_x, 'obj', stats.objval);
        end

        %% =================================================================
        % Recover duals by fixing y and solving a continuous QP
        function [result, stats] = solve_qp_fixed_y(obj, qpec, ~)
            % Build QP with all variables continuous, fix complementarity side by y_star,
            % and recover CasADi-style multipliers lam_g and reduced-costs lam_x.

            n  = qpec.dims.n_x;
            m  = qpec.dims.n_a;
            nc = qpec.dims.n_comp;

            % --- Continuous QP model ---
            M = struct();
            M.modelsense = 'min';
            M.vtype = repmat('C', n, 1);
            M.lb = qpec.lbw;
            M.ub = qpec.ubw;
            M.Q  = sparse(qpec.Q);
            M.obj = qpec.q;

            % --- Regular constraints with row map for lam_g reconstruction ---
            Arows = {}; rhs = {}; sense = {};
            map = struct('Aeq',nan(m,1),'Aub',nan(m,1),'Alb',nan(m,1));
            A = qpec.A; b = qpec.b; lba = qpec.lba; uba = qpec.uba;

            rcount = 0;
            for i = 1:m
                bi = b(i);
                if isfinite(lba(i)) && isfinite(uba(i)) && lba(i)==uba(i)
                    rcount = rcount + 1;
                    Arows{end+1,1} =  A(i,:); rhs{end+1,1} =  uba(i) - bi; sense{end+1,1} = '=';
                    map.Aeq(i) = rcount;
                else
                    if isfinite(uba(i))
                        rcount = rcount + 1;
                        Arows{end+1,1} =  A(i,:); rhs{end+1,1} =  uba(i) - bi; sense{end+1,1} = '<';
                        map.Aub(i) = rcount;
                    end
                    if isfinite(lba(i))
                        rcount = rcount + 1;
                        Arows{end+1,1} = -A(i,:); rhs{end+1,1} = -(lba(i) - bi); sense{end+1,1} = '<';
                        map.Alb(i) = rcount;
                    end
                end
            end

            % --- Complementarity fixed by y* ---
            y = round(qpec.y0(:));
            G = qpec.G; H = qpec.H; g = qpec.g; h = qpec.h;
            for i = 1:nc
                if y(i) == 0
                    % G_i(w) == -g_i ;  H_i(w) >= -h_i  -> -H_i w <= h_i
                    Arows{end+1,1} =  G(i,:); rhs{end+1,1} = -g(i); sense{end+1,1} = '=';
                    Arows{end+1,1} = -H(i,:); rhs{end+1,1} =  h(i); sense{end+1,1} = '<';
                else
                    % H_i(w) == -h_i ;  G_i(w) >= -g_i  -> -G_i w <= g_i
                    Arows{end+1,1} =  H(i,:); rhs{end+1,1} = -h(i); sense{end+1,1} = '=';
                    Arows{end+1,1} = -G(i,:); rhs{end+1,1} =  g(i); sense{end+1,1} = '<';
                end
            end

            M.A     = sparse(vertcat(Arows{:}));
            M.rhs   = vertcat(rhs{:});
            M.sense = char(vertcat(sense{:}));

            % optonally warm start qp?
            if obj.options.warmstart_qp
                M.start = qpec.w0(:);
            end

            % Ensure convex QP path returns duals
            M.Q = 0.5*(M.Q + M.Q'); % symmetrize
            params_qp = struct('OutputFlag',1,'Method',-1,'Crossover',1);

            % --- Solve QP ---
            t0   = tic;
            gr   = gurobi(M, params_qp);
            cpu  = toc(t0);
            stats = struct('status', obj.field_or(gr,'status','unknown'), 'cpu_time', cpu);

            lam_x = zeros(n,1);
            lam_g = zeros(m,1);
            lam_x = zeros(n,1);
            lam_G = zeros(nc,1);
            lam_H = zeros(nc,1);
            lam_G0 = zeros(nc,1);
            lam_H0 = zeros(nc,1);
            result = gr;
            stats.success = isfield(gr,'status') && (strcmpi(gr.status,'OPTIMAL') || strcmpi(gr.status,'SUBOPTIMAL'));
            %
            % if ~isfield(gr,'pi') || ~strcmp(stats.status,'OPTIMAL')
            %     return;
            % end
            % recover multipliers
            if stats.success
                pi = gr.pi(:);

                % --- lam_g reconstruction (CasADi style) ---
                for i = 1:m
                    if ~isnan(map.Aeq(i))
                        lam_g(i) = pi(map.Aeq(i));
                    else
                        ub_i = map.Aub(i); lb_i = map.Alb(i);
                        if ~isnan(ub_i) && ~isnan(lb_i)
                            lam_g(i) = pi(ub_i) - pi(lb_i);
                        elseif ~isnan(ub_i)
                            lam_g(i) = pi(ub_i);
                        elseif ~isnan(lb_i)
                            lam_g(i) = -pi(lb_i);
                        end
                    end
                end

                % --- lam_x from reduced costs ---
                if isfield(gr,'rc')
                    lam_x = gr.rc(:);
                end

                % mpec multiplieres:
                x = gr.x(:);
                g_eval = qpec.g + qpec.G*x;
                h_eval = qpec.h + qpec.H*x;

                tol = 1e-8;

                % classify activeness directly
                for i = 1:nc
                    if abs(g_eval(i)) <= tol  % active
                        lam_G(i) = max(0, - (g_eval(i))); % or just keep placeholder, depends on need
                    else
                        lam_G(i) = 0;
                    end
                    if abs(h_eval(i)) <= tol  % active
                        lam_H(i) = max(0, - (h_eval(i)));
                    else
                        lam_H(i) = 0;
                    end
                    if abs(g_eval(i)) + abs(h_eval(i)) <= 1e-6
                        lam_G0(i) = lam_G(i);
                        lam_H0(i) = lam_H(i);
                        if lam_H0(i) < 0 || lam_G0(i) < 0
                            stats.success = 0;
                        end
                    end
                end
            else
                result.x = zeros(n,1);
                % result.x = qpec.w0;
            end
            result.lam_x = lam_x;
            result.lam_g = lam_g;
            result.y = qpec.y0; % pass existing y0
            result.lam_G = lam_G;
            result.lam_H = lam_H;
            result.lam_G0 = lam_G0;
            result.lam_H0 = lam_H0;
            % --- One-line log  ---
            persistent qp_log_header_printed
            if isempty(qp_log_header_printed)
                fprintf('%-10s %-12s %12s %12s %12s %10s\n', ...
                    'Solve','Status','ObjVal','||x||_inf','Runtime','Nodes');
                qp_log_header_printed = true;
            end
            objval = obj.field_or(gr,'objval',NaN);
            % safe logging values
            if isfield(gr,'x') && ~isempty(gr.x)
                xnorm = norm(gr.x(1:n), inf);   % or norm(gr.x,inf) if you prefer
            else
                xnorm = NaN;
            end
            if isfield(gr,'objval')
                objval = gr.objval;
            else
                objval = NaN;
            end
            fprintf('%-10s %-12s %12.4g %12.4g %12.4g %10s\n', 'Piece QP', stats.status, objval, xnorm, cpu, '-');

        end

        %% =================================================================
        % Helpers
        function result = parse_result(~, grb_result, qpec)
            % Extract x, y, objective (primal only)
            n  = qpec.dims.n_x;
            nc = qpec.dims.n_comp;
            result = struct('x',[],'y',[],'obj',NaN,'lam_g',[],'lam_x',[]);
            if isfield(grb_result,'x')
                result.x = grb_result.x(1:n);
                if nc>0 && numel(grb_result.x) >= n+nc
                    result.y = grb_result.x(n+(1:nc));
                else
                    result.y = zeros(nc,1);
                end
            else
                result.x = zeros(n,1);
                result.y = zeros(nc,1);
            end
            if isfield(grb_result,'objval')
                result.obj = grb_result.objval;
            end
        end

        function out = field_or(~, s, f, def)
            out = def; if isfield(s,f), out = s.(f); end
        end
    end
end
