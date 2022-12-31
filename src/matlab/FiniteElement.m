classdef FiniteElement < NosnocFormulationObject
    properties
        % Index vectors
        ind_x
        ind_v
        ind_theta
        ind_lam
        ind_mu
        ind_alpha
        ind_lambda_n
        ind_lambda_p
        ind_h
        ind_nu_lift
        ind_elastic
        ind_boundary % index of bundary value lambda and mu, TODO is this even necessary?

        ind_g_comp_path
        
        ctrl_idx
        fe_idx
        
        model
        settings
        dims
        
        prev_fe
    end

    properties(Dependent, SetAccess=private, Hidden)
        x
        v
        theta
        lam
        mu
        alpha
        lambda_n
        lambda_p
        h

        elastic

        lambda
    end

    methods
        function obj = FiniteElement(prev_fe, settings, model,  dims, ctrl_idx, fe_idx)
            import casadi.*
            obj@NosnocFormulationObject();
            obj.ind_x = cell(dims.n_s, 1);
            obj.ind_v = cell(dims.n_s, 1);
            obj.ind_theta = cell(dims.n_s, dims.n_sys);
            obj.ind_lam = cell(dims.n_s, dims.n_sys);
            obj.ind_mu = cell(dims.n_s, dims.n_sys);
            obj.ind_alpha = cell(dims.n_s, dims.n_sys);
            obj.ind_lambda_n = cell(dims.n_s, dims.n_sys);
            obj.ind_lambda_p = cell(dims.n_s, dims.n_sys);
            obj.ind_h = [];
            obj.ind_nu_lift = [];
            obj.ind_elastic = cell(dims.n_s, 1);
            obj.ind_boundary = [];

            obj.ctrl_idx = ctrl_idx;
            obj.fe_idx = fe_idx;
            
            obj.settings = settings;
            obj.model = model;
            obj.dims = dims;

            obj.prev_fe = prev_fe;

            h = SX.sym(['h_' num2str(ctrl_idx) '_' num2str(fe_idx)]);
            h_ctrl_stage = model.T/dims.N_stages;
            h0 = h_ctrl_stage / dims.N_finite_elements;
            ubh = (1 + settings.gamma_h) * h0;
            lbh = (1 - settings.gamma_h) * h0;
            obj.addVariable(h, 'h', lbh, ubh, h0);

            % RK stage stuff
            for ii = 1:dims.n_s
                % state / state derivative variables
                if (settings.irk_representation == IrkRepresentation.differential ||...
                    settings.irk_representation == IrkRepresentation.differential_lift_x)
                    obj.addVariable(SX.sym(['V_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii)'], dims.n_x),...
                                    'v',...
                                    -inf * ones(dims.n_x),...
                                    inf * ones(dims.n_x),...
                                    zeros(dims.n_x),...
                                    ii);
                end
                if (settings.irk_representation == IrkRepresentation.integral ||...
                    settings.irk_representation == IrkRepresentation.differential_lift_x)
                    obj.addVariable(SX.sym(['X_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii)], dims.n_x),...
                                    'x',...
                                    -inf * ones(dims.n_x),...
                                    inf * ones(dims.n_x),...
                                    model.x0,...
                                    ii);
                end
                % algebraic variables
                if settings.pss_mode == PssMode.Stewart
                    % add thetas
                    for ij = 1:dims.n_sys
                        obj.addVariable(SX.sym(['theta_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)], dims.n_f_sys(ij)),...
                                        'theta',...
                                        zeros(dims.n_f_sys(ij)),...
                                        inf * ones(dims.n_f_sys(ij)),...
                                        settings.initial_theta * ones(dims.n_f_sys(ij)),...
                                        ii,...
                                        ij);
                    end
                    % add lambdas
                    for ij = 1:dims.n_sys
                        obj.addVariable(SX.sym(['lambda_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)], dims.n_f_sys(ij)),...
                                        'lam',...
                                        -inf * ones(dims.n_f_sys(ij)),...
                                        inf * ones(dims.n_f_sys(ij)),...
                                        settings.initial_lambda * ones(dims.n_f_sys(ij)),...
                                        ii,...
                                        ij);
                    end
                    % add mu
                    for ij = 1:dims.n_sys
                        obj.addVariable(SX.sym(['mu_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)], 1),...
                                        'mu',...
                                        -inf * ones(1),...
                                        inf * ones(1),...
                                        settings.initial_mu * ones(1),...
                                        ii,...
                                        ij);
                    end
                elseif settings.pss_mode == PssMode.Step
                    % add alpha
                    for ij = 1:dims.n_sys
                        obj.addVariable(SX.sym(['alpha_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)], dims.n_c_sys(ij)),...
                                        'alpha',...
                                        zeros(dims.n_c_sys(ij)),...
                                        ones(dims.n_c_sys(ij)),...
                                        settings.initial_theta * ones(dims.n_c_sys(ij)),...
                                        ii,...
                                        ij);
                    end
                    % add lambda_n and lambda_p
                    for ij = 1:dims.n_sys
                        obj.addVariable(SX.sym(['lambda_n_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)], dims.n_c_sys(ij)),...
                                        'lambda_n',...
                                        zeros(dims.n_c_sys(ij)),...
                                        inf * ones(dims.n_c_sys(ij)),...
                                        settings.initial_lambda_0 * ones(dims.n_c_sys(ij)),...
                                        ii,...
                                        ij);
                        obj.addVariable(SX.sym(['lambda_p_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)], dims.n_c_sys(ij)),...
                                        'lambda_p',...
                                        zeros(dims.n_c_sys(ij)),...
                                        inf * ones(dims.n_c_sys(ij)),...
                                        settings.initial_lambda_1 * ones(dims.n_c_sys(ij)),...
                                        ii,...
                                        ij);
                    end
                end
            end
            % Add right boundary points if needed
            if ~settings.right_boundary_point_explicit
                if settings.pss_mode == PssMode.Stewart
                    % add lambdas
                    for ij = 1:dims.n_sys
                        obj.addVariable(SX.sym(['lambda_' num2str(ctrl_idx) '_' num2str(fe_idx) '_end_' num2str(ij)], dims.n_f_sys(ij)),...
                                        'lam',...
                                        zeros(dims.n_f_sys(ij)),...
                                        inf * ones(dims.n_f_sys(ij)),...
                                        settings.initial_lambda * ones(dims.n_f_sys(ij)),...
                                        dims.n_s,...
                                        ij);
                    end
                    % add mu
                    for ij = 1:dims.n_sys
                        obj.addVariable(SX.sym(['mu_' num2str(ctrl_idx) '_' num2str(fe_idx) '_end_' num2str(ij)], 1),...
                                        'mu',...
                                        -inf * ones(1),...
                                        inf * ones(1),...
                                        settings.initial_mu * ones(1),...
                                        dims.n_s,...
                                        ij);
                    end
                elseif settings.pss_mode == PssMode.Step
                    % add lambda_n
                    for ij = 1:dims.n_sys
                        obj.addVariable(SX.sym(['lambda_n_' num2str(ctrl_idx) '_' num2str(fe_idx) '_end_' num2str(ij)], dims.n_c_sys(ij)),...
                                        'lambda_n',...
                                        zeros(dims.n_c_sys(ij)),...
                                        inf * ones(dims.n_c_sys(ij)),...
                                        settings.initial_lambda_0 * ones(dims.n_c_sys(ij)),...
                                        dims.n_s,...
                                        ij);
                    end
                    % add lambda_p
                    for ij = 1:dims.n_sys
                        obj.addVariable(SX.sym(['lambda_p_' num2str(ctrl_idx) '_' num2str(fe_idx) '_end_' num2str(ij)], dims.n_c_sys(ij)),...
                                        'lambda_p',...
                                        zeros(dims.n_c_sys(ij)),...
                                        inf * ones(dims.n_c_sys(ij)),...
                                        settings.initial_lambda_1 * ones(dims.n_c_sys(ij)),...
                                        dims.n_s,...
                                        ij);
                    end
                end
            end
            if (~settings.right_boundary_point_explicit ||...
                settings.irk_representation == IrkRepresentation.differential)
                % add final X variables
                obj.addVariable(SX.sym(['X_end_' num2str(ctrl_idx) '_' num2str(fe_idx)], dims.n_x),...
                                'x',...
                                -inf * ones(dims.n_x),...
                                inf * ones(dims.n_x),...
                                model.x0,...
                                ii+1);
            end

            if strcmpi(settings.step_equilibration, 'direct_homotopy_lift')
                obj.addVariable(SX.sym(['nu_lift_' num2str(ctrl_idx) '_' num2str(fe_idx)], 1),...
                                'nu_lift',...
                                1,...
                                -inf,...
                                inf);
            end
        end
        
        function lambda = get.lambda(obj)
            grab = @(l, ln, lp) vertcat(obj.w(l), obj.w(ln), obj.w(lp));

            lambda = cellfun(grab, obj.ind_lam, obj.ind_lambda_p, obj.ind_lambda_n, 'UniformOutput', false);
        end

        function theta = get.theta(obj)
            grab = @(t, a) vertcat(obj.w(t), obj.w(a), ones(size(a)) - obj.w(a));

            theta = cellfun(grab, obj.ind_theta, obj.ind_alpha, 'UniformOutput', false);
        end

        function h = get.h(obj)
            h = obj.w(obj.ind_h);
        end

        function x = get.x(obj)
            x = cellfun(@(x) obj.w(x), obj.ind_x, 'UniformOutput', false);
        end

        function v = get.v(obj)
            v = cellfun(@(v) obj.w(v), obj.ind_v, 'UniformOutput', false);
        end

        function elastic = get.elastic(obj)
            elastic = cellfun(@(ell) obj.w(ell), obj.ind_elastic, 'UniformOutput', false);
        end
        
        function sum_lambda = sumLambda(obj, varargin)
            import casadi.*
            p = inputParser();
            p.FunctionName = 'sumLambda';
            
            % TODO: add checks.
            addRequired(p, 'obj');
            addOptional(p, 'sys',[]);
            parse(p, obj, varargin{:});

            if ismember('sys', p.UsingDefaults)
                lambda = obj.lambda;
                lambdas = arrayfun(@(row) vertcat(lambda{row, :}), 1:size(lambda,1), 'UniformOutput', false);
                lambdas = [lambdas, {vertcat(obj.prev_fe.lambda{end,:})}];
            else
                lambdas = obj.lambda(:,p.Results.sys).';
                lambdas = [lambdas, obj.prev_fe.lambda(end,p.Results.sys)];
            end
            sum_lambda = sum([lambdas{:}], 2);
        end

        function sum_theta = sumTheta(obj, varargin)
            p = inputParser();
            p.FunctionName = 'sumLambda';
            
            % TODO: add checks.
            addRequired(p, 'obj');
            addOptional(p, 'sys',[]);
            parse(p, obj, varargin{:});

            if ismember('sys', p.UsingDefaults)
                theta = obj.theta;
                thetas = arrayfun(@(row) vertcat(theta{row, :}), 1:size(theta,1), 'UniformOutput', false);
            else
                thetas = obj.theta{:,p.Results.sys}.';
            end
            sum_theta = sum([thetas{:}], 2);
        end

        function sum_elastic = sumElastic(obj)
            elastic = obj.elastic;
            sum_elastic = sum([elastic{:}]);
        end

        function z = rkStageZ(obj, stage)
            import casadi.*

            idx = [[obj.ind_theta{stage, :}],...
                   [obj.ind_lam{stage, :}],...
                   [obj.ind_mu{stage, :}],...
                   [obj.ind_alpha{stage, :}],...
                   [obj.ind_lambda_n{stage, :}],...
                   [obj.ind_lambda_p{stage, :}]];

            z = obj.w(idx);
        end

        function forwardSimulation(obj, ocp, Uk, s_sot_k)
            model = obj.model;
            settings = obj.settings;
            dims = obj.dims;

            if settings.irk_representation == IrkRepresentation.integral
                X_ki = obj.x;
                Xk_end = settings.D_irk(1) * obj.prev_fe.x{end};
            elseif settings.irk_representation == IrkRepresentation.differential
                X_ki = {};
                for j = 1:dims.n_s
                    x_temp = obj.prev_fe.x{end};
                    for r = 1:dims.n_s
                        x_temp = x_temp + obj.h*settings.A_irk(j,r)*obj.v{r};
                    end
                    X_ki = [X_ki{:} x_temp];
                end
                X_ki = [X_ki{:}, obj.x{end}];
                Xk_end = obj.prev_fe.x{end};
            elseif settings.irk_representation == IrkRepresentation.differential_lift_x
                X_ki = {};
                for j = 1:dims.n_s
                    x_temp = obj.prev_fe.x{end};
                    for r = 1:dims.n_s
                        x_temp = x_temp + obj.h*settings.A_irk(j,r)*obj.v{r};
                    end
                    X_ki = [X_ki{:} x_temp];
                    obj.addConstraint(self.x{j}-x_temp)
                end
                X_ki = [X_ki{:}, obj.x{end}];
                Xk_end = obj.prev_fe.x{end};
            end

            for j = 1:dims.n_s
                % Multiply by s_sot_k which is 1 if not using speed of time variables
                [fj, qj] = model.f_x_fun(X_ki{j}, obj.rkStageZ(j), Uk);
                fj = s_sot_k*fj;
                qj = s_sot_k*qj;
                gj = s_sot_k*model.g_z_all_fun(X_ki{j}, obj.rkStageZ(j), Uk);
                
                obj.addConstraint(gj);
                if settings.irk_representation == IrkRepresentation.integral
                    xj = settings.C_irk(1, j+1) * obj.prev_fe.x{end};
                    for r = 1:dims.n_s
                        xj = xj + settings.C_irk(r+1, j+1) * X_ki{r};
                    end
                    Xk_end = Xk_end + settings.D_irk(j+1) * X_ki{j};
                    obj.addConstraint(obj.h * fj - xj);
                    obj.cost = obj.cost + settings.B_irk(j+1) * obj.h * qj;
                else
                    Xk_end = Xk_end + obj.h * settings.b_irk(j) * obj.v{j};
                    obj.add_constraint(fj - obj.v{j})
                    obj.cost = obj.cost + settings.b_irk(j) * obj.h * qj
                end
            end

            % nonlinear inequality.
            % TODO: do this cleaner
            for j=1:dims.n_s-settings.right_boundary_point_explicit
                % TODO: there has to be a better way to do this.
                if settings.g_ineq_constraint && settings.g_ineq_at_stg
                    obj.addContsraint(model.g_ineq_fun(X_ki{j},Uk), model.g_ineq_lb, model.g_ineq_ub);
                end
            end
            if (settings.g_ineq_constraint &&...
                (obj.fe_idx == dims.N_finite_elements(obj.ctrl_idx) || settings.g_ineq_at_fe))
                obj.addContsraint(model.g_ineq_fun(X_ki{end},Uk), model.g_ineq_lb, model.g_ineq_ub);
            end

            % path complementarities
            % TODO: do this cleaner
            for j=1:dims.n_s-settings.right_boundary_point_explicit
                % TODO: there has to be a better way to do this.
                if settings.g_comp_path_constraint && settings.g_ineq_at_stg
                    obj.addContsraint(model.g_comp_path_fun(X_ki{j},Uk), model.g_ineq_lb, model.g_ineq_ub);
                end
            end
            if (settings.g_comp_path_constraint &&...
                (obj.fe_idx == dims.N_finite_elements(obj.ctrl_idx) || settings.g_ineq_at_fe))
                obj.addContsraint(model.g_ineq_fun(X_ki{end},Uk), model.g_ineq_lb, model.g_ineq_ub);
            end

            % end constraints
            if (~settings.right_boundary_point_explicit ||...
                settings.irk_representation == IrkRepresentation.differential)
                self.addConstraint(Xk_end - obj.x{end});
            end
            if (~settings.right_boundary_point_explicit &&...
                settings.use_fesd &&...
                obj.fe_idx < dims.N_finite_elements(obj.ctrl_idx)) % TODO make this handle different numbers of FE
                
                % TODO verify this.
                obj.addConstraint(model.g_z_switching_fun(obj.x{end}, self.rkStageZ(dims.n_s+1), Uk));
            end
        end

        function createComplementarityConstraints(obj, sigma_p, s_elastic)
            import casadi.*           
            model = obj.model;
            settings = obj.settings;
            dims = obj.dims;

            g_cross_comp = SX([]);
            % TODO Implement other modes
            if ~settings.use_fesd
                % TODO vectorize?
                for j=1:dims.n_s
                    theta = vertcat(obj.theta{j,:});
                    lambda = vertcat(obj.lambda{j,:});
                    g_cross_comp = vertcat(g_cross_comp, dot(theta,lambda));
                end
            elseif settings.cross_comp_mode == 1 
                % TODO vectorize?
                theta = obj.theta;
                lambda = obj.lambda;
                % Complement within FE
                for j=1:dims.n_s
                    for jj = 1:dims.n_s
                        for r=1:dims.n_sys 
                            g_cross_comp = vertcat(g_cross_comp, diag(theta{j,r})*lambda{jj,r});
                        end
                    end
                end
                for j=1:dims.n_s
                    for r=1:dims.n_sys
                        g_cross_comp = vertcat(g_cross_comp, diag(theta{j,r})*obj.prev_fe.lambda{end,r});
                    end
                end
            elseif settings.cross_comp_mode == 2
                theta = obj.theta;
                lambda = obj.lambda;
                for r=1:dims.n_sys
                    for j=1:dims.n_s
                        for jj=1:dims.n_s
                            g_cross_comp = vertcat(g_cross_comp, dot(theta{j,r},lambda{jj,r}));
                        end
                    end
                end
                for r=1:dims.n_sys
                    for j=1:dims.n_sys
                        g_cross_comp = vertcat(g_cross_comp, dot(theta{j,r}, obj.prev_fe.lambda{end,r}));
                    end
                end
            elseif settings.cross_comp_mode == 3
                % TODO vectorize?
                theta = obj.theta;
                for j=1:dims.n_s
                    for r=1:dims.n_sys 
                        sum_lambda = obj.sumLambda(r);
                        g_cross_comp = vertcat(g_cross_comp, diag(theta{j,r})*sum_lambda);
                    end
                end
            elseif settings.cross_comp_mode == 4
                lambda = obj.lambda;
                for r=1:dims.n_sys
                    sum_theta = obj.sumTheta(r);
                    for j=1:dims.n_s
                        g_cross_comp = vertcat(g_cross_comp, diag(sum_theta)*lambda{j,r});
                    end
                    g_cross_comp = vertcat(g_cross_comp, diag(sum_theta)*obj.prev_fe.lambda{end,r});
                end

            elseif settings.cross_comp_mode == 5
                theta = obj.theta;
                for j=1:dims.n_s
                    for r=1:dims.n_sys
                        sum_lambda = obj.sumLambda(r);
                        g_cross_comp = vertcat(g_cross_comp, dot(theta{j,r}, sum_lambda));
                    end
                end
            elseif settings.cross_comp_mode == 6
                for r=1:dims.n_sys
                    sum_theta = obj.sumTheta(r);
                    for j=1:dims.n_s
                        g_cross_comp = vertcat(g_cross_comp, dot(sum_theta,lambda{j, r}));
                    end
                    g_cross_comp = vertcat(g_cross_comp, dot(sum_theta,obj.prev_fe.lambda{end,r}));
                end
            elseif settings.cross_comp_mode == 7
                error('TODO: not implemented');
            elseif settings.cross_comp_mode == 8
                error('TODO: not implemented');
            elseif settings.cross_comp_mode == 9
                error('TODO: not implemented');
            elseif settings.cross_comp_mode == 10
                error('TODO: not implemented');
            elseif settings.cross_comp_mode == 11
                error('TODO: not implemented');
            elseif settings.cross_comp_mode == 12
                error('TODO: not implemented');
            end

           
            n_cross_comp = length(g_cross_comp);
            %
            if ismember(settings.mpcc_mode, MpccMode.elastic_ell_1)
                for ii=1:n_s
                    obj.addVariable(SX.sym(['s_elastic_' num2str(obj.ctrl_idx) '_' num2str(obj.fe_idx) '_' num2str(ii), ], n_cross_comp),...
                                    'elastic',...
                                    s_elastic_0*ones(n_all_comp_j,1),...
                                    s_elastic_min*ones(n_all_comp_j,1),...
                                    s_elastic_max*ones(n_all_comp_j,1),...
                                    ii)
                end
            end
            
            % Do MPCC formulation
            % TODO this should be done on the problem level, and take into account passed in vars. fine for now.
            if settings.mpcc_mode == 'Scholtes_ineq'
                g_cross_comp = g_cross_comp - sigma_p;
                g_cross_comp_ub = zeros(n_cross_comp,1);
                g_cross_comp_lb = -inf * ones(n_cross_comp,1);
            elseif settings.mpcc_mode == 'Scholtes_eq'
                g_cross_comp = g_cross_comp - sigma_p;
                g_cross_comp_ub = zeros(n_cross_comp,1);
                g_cross_comp_lb = zeros(n_cross_comp,1);
            elseif settings.mpcc_mode == 'ell_1_penalty'
                if settings.objective_scaling_direct
                    obj.cost = obj.cost + (1/sigma_p)*sum(g_cross_comp);
                else
                    obj.cost = sigma_p*obj.cost + sum(g_cross_comp);
                end
            elseif settings.mpcc_mode == 'elastic_ineq'
                g_cross_comp = g_cross_comp - s_elastic;
                g_cross_comp_ub = zeros(n_cross_comp,1);
                g_cross_comp_lb = -inf * ones(n_cross_comp,1);
            elseif settings.mpcc_mode == 'elastic_eq'
                g_cross_comp = g_cross_comp - s_elastic;
                g_cross_comp_ub = zeros(n_cross_comp,1);
                g_cross_comp_lb = zeros(n_cross_comp,1);
            elseif settings.mpcc_mode == 'elastic_two_sided'
                g_cross_comp = g_cross_comp - s_elastic;
                g_cross_comp_ub = zeros(n_cross_comp,1);
                g_cross_comp_lb = -inf * ones(n_cross_comp,1);

                g_cross_comp = [g_cross_comp;g_all_comp_j+s_elastic*ones(n_cross_comp,1)];
                g_cross_comp_ub = [g_cross_comp_ub; inf*ones(n_cross_comp,1)];
                g_cross_comp_lb = [g_cross_comp_lb;  zeros(n_cross_comp,1)];
                % TODO these should be merged probably
            elseif settings.mpcc_mode == 'elastic_ell_1_ineq'
                g_cross_comp = g_cross_comp - s_elastic;
                g_cross_comp_ub = zeros(n_cross_comp,1);
                g_cross_comp_lb = -inf * ones(n_cross_comp,1);
            elseif settings.mpcc_mode == 'elastic_ell_1_eq'
                g_cross_comp = g_cross_comp - s_elastic;
                g_cross_comp_ub = zeros(n_cross_comp,1);
                g_cross_comp_lb = zeros(n_cross_comp,1);
            elseif settings.mpcc_mode == 'elastic_ell_1_two_sided'
                g_cross_comp = g_cross_comp - s_elastic;
                g_cross_comp_ub = zeros(n_cross_comp,1);
                g_cross_comp_lb = -inf * ones(n_cross_comp,1);

                g_cross_comp = [g_cross_comp;g_all_comp_j+s_elastic*ones(n_cross_comp,1)];
                g_cross_comp_ub = [g_cross_comp_ub; inf*ones(n_cross_comp,1)];
                g_cross_comp_lb = [g_cross_comp_lb;  zeros(n_cross_comp,1)];
            end

            if settings.mpcc_mode ~= 'ell_1_penalty'
                obj.addConstraint(g_cross_comp, g_cross_comp_lb, g_cross_comp_ub);
            end
        end

        function stepEquilibration(obj)
            import casadi.*
            model = obj.model;
            settings = obj.settings;
            dims = obj.dims;

            % TODO implement other modes!
            if settings.use_fesd && obj.fe_idx > 1
                delta_h_ki = obj.h - obj.prev_fe.h;
            end
            if settings.step_equilibration == StepEquilibrationMode.heuristic_mean
                h_fe = model.T / (dims.N_stages * dims.N_finite_elements);
                obj.cost = obj.cost + settings.rho_h * (obj.h - h_fe).^2;
            elseif settings.step_equilibration == StepEquilibrationMode.heuristic_delta
                obj.cost = obj.cost + settings.rho_h * delta_h_ki.^2;
            elseif settings.step_equilibration == StepEquilibrationMode.l2_relaxed_scaled
                eta_k = obj.prev_fe.sumLambda() * obj.sumLambda() + obj.prev_fe.sumTheta() * obj.sumTheta();
                nu_k = 1;
                for jjj=1:length(eta_k)
                    nu_k = nu_k * eta_k(jjj)
                end
                obj.cost = obj.cost + settings.rho_h * tanh(nu_k/settings.step_equilibration_sigma) * delta_h_ki.^2;
            elseif settings.step_equilibration == StepEquilibrationMode.l2_relaxed
                eta_k = obj.prev_fe.sumLambda() * obj.sumLambda() + obj.prev_fe.sumTheta() * obj.sumTheta();
                nu_k = 1;
                for jjj=1:length(eta_k)
                    nu_k = nu_k * eta_k(jjj)
                end
                obj.cost = obj.cost + settings.rho_h * nu_k * delta_h_ki.^2
            else
                error("Step equilibration mode not implemented");
            end
        end
    end
end

