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
        ind_beta
        ind_theta_step
        ind_z

        ctrl_idx
        fe_idx

        model
        settings
        dims

        u

        T_final

        prev_fe
    end

    properties(Dependent, SetAccess=private, Hidden)
        x
        v
        lam
        mu
        alpha
        lambda_n
        lambda_p
        nu_lift
        h

        elastic
        
        nu_vector
    end

    properties(SetAccess=private)
        lambda
        theta
    end

    methods
        function obj = FiniteElement(prev_fe, settings, model, dims, ctrl_idx, fe_idx, varargin)
            import casadi.*
            obj@NosnocFormulationObject();

            p = inputParser();
            p.FunctionName = 'FiniteElement';

            % TODO: add checks.
            addRequired(p, 'prev_fe');
            addRequired(p, 'settings');
            addRequired(p, 'model');
            addRequired(p, 'dims');
            addRequired(p, 'ctrl_idx');
            addRequired(p, 'fe_idx');
            addOptional(p, 'T_final',[]);
            parse(p, prev_fe, settings, model, dims, ctrl_idx, fe_idx, varargin{:});

            if settings.right_boundary_point_explicit
                rbp_allowance = 0;
            else
                rbp_allowance = 1;
            end

            obj.ind_x = cell(dims.n_s+rbp_allowance, 1);
            obj.ind_v = cell(dims.n_s, 1);
            obj.ind_theta = cell(dims.n_s+rbp_allowance, dims.n_sys);
            obj.ind_lam = cell(dims.n_s+rbp_allowance, dims.n_sys);
            obj.ind_mu = cell(dims.n_s+rbp_allowance, dims.n_sys);
            obj.ind_alpha = cell(dims.n_s+rbp_allowance, dims.n_sys);
            obj.ind_lambda_n = cell(dims.n_s+rbp_allowance, dims.n_sys);
            obj.ind_lambda_p = cell(dims.n_s+rbp_allowance, dims.n_sys);
            obj.ind_theta_step = cell(dims.n_s+rbp_allowance, 1);
            obj.ind_beta = cell(dims.n_s+rbp_allowance, 1);
            obj.ind_z = cell(dims.n_s+rbp_allowance, 1);
            obj.ind_h = [];
            obj.ind_elastic = [];
            obj.ind_boundary = [];

            obj.ctrl_idx = ctrl_idx;
            obj.fe_idx = fe_idx;

            obj.settings = settings;
            obj.model = model;
            obj.dims = dims;

            obj.prev_fe = prev_fe;

            if settings.use_fesd
                h = define_casadi_symbolic(settings.casadi_symbolic_mode, ['h_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1)]);
                h_ctrl_stage = model.T/dims.N_stages;
                h0 = h_ctrl_stage / dims.N_finite_elements(ctrl_idx);
                ubh = (1 + settings.gamma_h) * h0;
                lbh = (1 - settings.gamma_h) * h0;
                if settings.time_rescaling && ~settings.use_speed_of_time_variables
                    % if only time_rescaling is true, speed of time and step size all lumped together, e.g., \hat{h}_{k,i} = s_n * h_{k,i}, hence the bounds need to be extended.
                    ubh = (1+settings.gamma_h)*h0*settings.s_sot_max;
                    lbh = (1-settings.gamma_h)*h0/settings.s_sot_min;
                end
                obj.addVariable(h, 'h', lbh, ubh, h0);
            end
            obj.T_final = p.Results.T_final;

            % RK stage stuff
            for ii = 1:dims.n_s
                % state / state derivative variables
                if (settings.irk_representation == IrkRepresentation.differential ||...
                        settings.irk_representation == IrkRepresentation.differential_lift_x)
                    v = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                        ['V_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)'],...
                        dims.n_x);
                    obj.addVariable(v,...
                        'v',...
                        -inf * ones(dims.n_x,1),...
                        inf * ones(dims.n_x,1),...
                        zeros(dims.n_x,1),...
                        ii);
                end
                if (settings.irk_representation == IrkRepresentation.integral ||...
                        settings.irk_representation == IrkRepresentation.differential_lift_x)
                    if settings.x_box_at_stg
                        lbx = model.lbx;
                        ubx = model.ubx;
                    elseif settings.x_box_at_fe && ii == dims.n_s && settings.right_boundary_point_explicit
                        lbx = model.lbx;
                        ubx = model.ubx;
                    elseif fe_idx == dims.N_finite_elements(ctrl_idx) && ii == dims.n_s && settings.right_boundary_point_explicit
                        lbx = model.lbx;
                        ubx = model.ubx;
                    else
                        lbx = -inf * ones(dims.n_x, 1);
                        ubx = inf * ones(dims.n_x, 1);
                    end
                    x = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                        ['X_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                        dims.n_x);
                    obj.addVariable(x,...
                        'x',...
                        lbx,...
                        ubx,...
                        model.x0,...
                        ii);
                end
                % algebraic variables
                if settings.dcs_mode == DcsMode.Stewart
                    % add thetas
                    for ij = 1:dims.n_sys
                        theta = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['theta_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            dims.n_f_sys(ij));
                        obj.addVariable(theta,...
                            'theta',...
                            zeros(dims.n_f_sys(ij), 1),...
                            inf * ones(dims.n_f_sys(ij), 1),...
                            settings.initial_theta * ones(dims.n_f_sys(ij), 1),...
                            ii,...
                            ij);
                    end
                    % add lambdas
                    for ij = 1:dims.n_sys
                        lam = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['lambda_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            dims.n_f_sys(ij));
                        obj.addVariable(lam,...
                            'lam',...
                            zeros(dims.n_f_sys(ij), 1),...
                            inf * ones(dims.n_f_sys(ij), 1),...
                            settings.initial_lambda * ones(dims.n_f_sys(ij), 1),...
                            ii,...
                            ij);
                    end
                    % add mu
                    for ij = 1:dims.n_sys
                        mu = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['mu_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            1);
                        obj.addVariable(mu,...
                            'mu',...
                            -inf,...
                            inf,...
                            settings.initial_mu,...
                            ii,...
                            ij);
                    end
                elseif settings.dcs_mode == DcsMode.Step
                    % add alpha
                    for ij = 1:dims.n_sys
                        alpha = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['alpha_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            dims.n_c_sys(ij));
                        obj.addVariable(alpha,...
                            'alpha',...
                            zeros(dims.n_c_sys(ij), 1),...
                            ones(dims.n_c_sys(ij), 1),...
                            settings.initial_theta * ones(dims.n_c_sys(ij), 1),...
                            ii,...
                            ij);
                    end
                    % add lambda_n and lambda_p
                    for ij = 1:dims.n_sys
                        lambda_n = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['lambda_n_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            dims.n_c_sys(ij));
                        lambda_p = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['lambda_p_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            dims.n_c_sys(ij));
                        obj.addVariable(lambda_n,...
                            'lambda_n',...
                            zeros(dims.n_c_sys(ij), 1),...
                            inf * ones(dims.n_c_sys(ij), 1),...
                            settings.initial_lambda_0 * ones(dims.n_c_sys(ij), 1),...
                            ii,...
                            ij);
                        obj.addVariable(lambda_p,...
                            'lambda_p',...
                            zeros(dims.n_c_sys(ij), 1),...
                            inf * ones(dims.n_c_sys(ij), 1),...
                            settings.initial_lambda_1 * ones(dims.n_c_sys(ij), 1),...
                            ii,...
                            ij);
                    end
                    % TODO: Clean this up (maybe as a function to reduce the indent level.
                    if settings.time_freezing_inelastic
                        %                         if ~settings.pss_lift_step_functions
                        theta_step = define_casadi_symbolic(settings.casadi_symbolic_mode, ['theta_step_'  num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],dims.n_theta_step);
                        beta = define_casadi_symbolic(settings.casadi_symbolic_mode, ['beta_'  num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],dims.n_beta);
                        beta_guess = full(model.g_lift_beta_fun(settings.initial_alpha*ones(dims.n_alpha,1)));
                        theta_step_guess = full(model.g_lift_theta_step_fun(settings.initial_alpha*ones(dims.n_alpha,1),beta_guess));
                        obj.addVariable(theta_step,...
                            'theta_step',...
                            -inf*ones(dims.n_theta_step,1),...
                            inf*ones(dims.n_theta_step,1),...
                            theta_step_guess,...
                            ii);
                        if settings.pss_lift_step_functions
                            obj.addVariable(beta,...
                                'beta',...
                                -inf*ones(dims.n_beta,1),...
                                inf*ones(dims.n_beta,1),...
                                beta_guess,...
                                ii);
                        end
                    end
                end

                % add user variables
                z = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                    ['z_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                    dims.n_z);
                obj.addVariable(z,...
                    'z',...
                    model.lbz,...
                    model.ubz,...
                    model.z0,...
                    ii);
            end
            % Add right boundary points if needed
            if ~settings.right_boundary_point_explicit
                if settings.dcs_mode == DcsMode.Stewart
                    % add lambdas
                    for ij = 1:dims.n_sys
                        lam = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['lambda_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_end_' num2str(ij)],...
                            dims.n_f_sys(ij));
                        obj.addVariable(lam,...
                            'lam',...
                            zeros(dims.n_f_sys(ij), 1),...
                            inf * ones(dims.n_f_sys(ij), 1),...
                            settings.initial_lambda * ones(dims.n_f_sys(ij), 1),...
                            dims.n_s+1,...
                            ij);
                    end

                    % add mu
                    for ij = 1:dims.n_sys
                        mu = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['mu_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_end_' num2str(ij)],...
                            1);
                        obj.addVariable(mu,...
                            'mu',...
                            -inf * ones(1),...
                            inf * ones(1),...
                            settings.initial_mu * ones(1),...
                            dims.n_s+1,...
                            ij);
                    end

                elseif settings.dcs_mode == DcsMode.Step
                    % add lambda_n
                    for ij = 1:dims.n_sys
                        lambda_n = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['lambda_n_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_end_' num2str(ij)],...
                            dims.n_c_sys(ij));
                        obj.addVariable(lambda_n,...
                            'lambda_n',...
                            zeros(dims.n_c_sys(ij), 1),...
                            inf * ones(dims.n_c_sys(ij), 1),...
                            settings.initial_lambda_0 * ones(dims.n_c_sys(ij), 1),...
                            dims.n_s+1,...
                            ij);
                    end

                    % add lambda_p
                    for ij = 1:dims.n_sys
                        lambda_p = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                            ['lambda_p_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_end_' num2str(ij)],...
                            dims.n_c_sys(ij));
                        obj.addVariable(lambda_p,...
                            'lambda_p',...
                            zeros(dims.n_c_sys(ij), 1),...
                            inf * ones(dims.n_c_sys(ij), 1),...
                            settings.initial_lambda_1 * ones(dims.n_c_sys(ij), 1),...
                            dims.n_s+1,...
                            ij);
                    end
                end
            end

            if (~settings.right_boundary_point_explicit ||...
                    settings.irk_representation == IrkRepresentation.differential)
                if settings.x_box_at_stg || settings.x_box_at_fe || fe_idx == dims.N_finite_elements(ctrl_idx)
                    lbx = model.lbx;
                    ubx = model.ubx;
                else
                    lbx = -inf * ones(dims.n_x);
                    ubx = inf * ones(dims.n_x);
                end

                % add final X variables
                x = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                    ['X_end_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1)],...
                    dims.n_x);
                obj.addVariable(x,...
                    'x',...
                    lbx,...
                    ubx,...
                    model.x0,...
                    dims.n_s+rbp_allowance);
            end
            if strcmpi(settings.step_equilibration, 'direct_homotopy_lift')
                nu_lift = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                    ['nu_lift_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1)],...
                    1);
                obj.addVariable(nu_lift,...
                    'nu_lift',...
                    1,...
                    -inf,...
                    inf);
            end

            % calculate lambda and theta
            grab = @(l, ln, lp) vertcat(obj.w(l), obj.w(ln), obj.w(lp));
            obj.lambda = cellfun(grab, obj.ind_lam, obj.ind_lambda_n, obj.ind_lambda_p, 'UniformOutput', false);
            grab = @(t, a) vertcat(obj.w(t), obj.w(a), ones(size(a))' - obj.w(a));
            obj.theta = cellfun(grab, obj.ind_theta, obj.ind_alpha, 'UniformOutput', false);
        end

        function h = get.h(obj)
            if obj.settings.use_fesd
                h = obj.w(obj.ind_h);
            elseif obj.settings.time_optimal_problem && ~obj.settings.use_speed_of_time_variables
                h = obj.T_final/(obj.dims.N_stages*obj.dims.N_finite_elements(obj.ctrl_idx));
            else
                h = obj.model.T/(obj.dims.N_stages*obj.dims.N_finite_elements(obj.ctrl_idx));
            end
        end

        function x = get.x(obj)
            x = cellfun(@(x) obj.w(x), obj.ind_x, 'UniformOutput', false);
        end

        function v = get.v(obj)
            v = cellfun(@(v) obj.w(v), obj.ind_v, 'UniformOutput', false);
        end

        function elastic = get.elastic(obj)
            elastic = obj.w(obj.ind_elastic);
        end

        function nu_lift = get.nu_lift(obj)
            nu_lift = obj.w(obj.ind_nu_lift);
        end

        function nu_vector = get.nu_vector(obj)
            import casadi.*
            if obj.settings.use_fesd && obj.fe_idx > 1
                eta_k = obj.prev_fe.sumLambda().*obj.sumLambda() + obj.prev_fe.sumTheta().*obj.sumTheta();
                nu_vector = 1;
                for jjj=1:length(eta_k)
                    nu_vector = nu_vector * eta_k(jjj);
                end
            else
                nu_vector = [];
            end
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
            p.FunctionName = 'sumTheta';

            % TODO: add checks.
            addRequired(p, 'obj');
            addOptional(p, 'sys',[]);
            parse(p, obj, varargin{:});
            %obj.theta
            if ismember('sys', p.UsingDefaults)
                theta = obj.theta;
                thetas = arrayfun(@(row) vertcat(theta{row, :}), 1:size(theta,1), 'UniformOutput', false);
            else
                thetas = obj.theta(:,p.Results.sys).';
            end
            sum_theta = sum([thetas{:}], 2);
        end

        function sum_elastic = sumElastic(obj)
            elastic = obj.elastic;
            sum_elastic = sum(elastic);
        end

        function z = rkStageZ(obj, stage)
            import casadi.*

            % TODO: theta_step/beta
            idx = [[obj.ind_theta{stage, :}],...
                [obj.ind_lam{stage, :}],...
                [obj.ind_mu{stage, :}],...
                [obj.ind_alpha{stage, :}],...
                [obj.ind_lambda_n{stage, :}],...
                [obj.ind_lambda_p{stage, :}],...
                [obj.ind_beta{stage}],...
                [obj.ind_theta_step{stage}],...
                [obj.ind_z{stage}]];

            z = obj.w(idx);
        end

        function forwardSimulation(obj, ocp, Uk, s_sot, p_stage)
            model = obj.model;
            settings = obj.settings;
            dims = obj.dims;

            obj.u = Uk;

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
                    X_ki = [X_ki {x_temp}];
                end
                X_ki = [X_ki, {obj.x{end}}];
                Xk_end = obj.prev_fe.x{end};
            elseif settings.irk_representation == IrkRepresentation.differential_lift_x
                X_ki = obj.x;
                Xk_end = obj.prev_fe.x{end};
                for j = 1:dims.n_s
                    x_temp = obj.prev_fe.x{end};
                    for r = 1:dims.n_s
                        x_temp = x_temp + obj.h*settings.A_irk(j,r)*obj.v{r};
                    end
                    obj.addConstraint(obj.x{j}-x_temp);
                end
            end

            for j = 1:dims.n_s
                % Multiply by s_sot_k which is 1 if not using speed of time variable
                [fj, qj] = model.f_x_fun(X_ki{j}, obj.rkStageZ(j), Uk, p_stage, model.v_global);
                fj = s_sot*fj;
                qj = s_sot*qj;
                gj = model.g_z_all_fun(X_ki{j}, obj.rkStageZ(j), Uk, p_stage, model.v_global);

                obj.addConstraint(gj);
                if settings.irk_representation == IrkRepresentation.integral
                    xj = settings.C_irk(1, j+1) * obj.prev_fe.x{end};
                    for r = 1:dims.n_s
                        xj = xj + settings.C_irk(r+1, j+1) * X_ki{r};
                    end
                    Xk_end = Xk_end + settings.D_irk(j+1) * X_ki{j};
                    obj.addConstraint(obj.h * fj - xj);
                    obj.cost = obj.cost + settings.B_irk(j+1) * obj.h * qj;
                    obj.objective = obj.objective + settings.B_irk(j+1) * obj.h * qj;
                else
                    Xk_end = Xk_end + obj.h * settings.b_irk(j) * obj.v{j};
                    obj.addConstraint(fj - obj.v{j});
                    obj.cost = obj.cost + settings.b_irk(j) * obj.h * qj;
                    obj.objective = obj.objective + settings.b_irk(j) * obj.h * qj;
                end
            end

            % nonlinear inequality.
            % TODO: do this cleaner
            if (model.g_path_constraint &&...
                    (obj.fe_idx == dims.N_finite_elements(obj.ctrl_idx) || settings.g_path_at_fe))
                obj.addConstraint(model.g_path_fun(obj.prev_fe.x{end},Uk,p_stage,model.v_global), model.g_path_lb, model.g_path_ub);
            end
            for j=1:dims.n_s-settings.right_boundary_point_explicit
                % TODO: there has to be a better way to do this.
                if model.g_path_constraint && settings.g_path_at_stg
                    obj.addConstraint(model.g_path_fun(X_ki{j},Uk,p_stage,model.v_global), model.g_path_lb, model.g_path_ub);
                end
            end

            % end constraints
            if (~settings.right_boundary_point_explicit ||...
                    settings.irk_representation == IrkRepresentation.differential)
                obj.addConstraint(Xk_end - obj.x{end});
            end
            if (~settings.right_boundary_point_explicit &&...
                    settings.use_fesd &&...
                    obj.fe_idx < dims.N_finite_elements(obj.ctrl_idx))

                % TODO verify this.
                obj.addConstraint(model.g_switching_fun(obj.x{end}, obj.rkStageZ(dims.n_s+1), Uk, p_stage));
            end
        end

        function createComplementarityConstraints(obj, sigma_p, s_elastic, p_stage)
            import casadi.*
            model = obj.model;
            settings = obj.settings;
            dims = obj.dims;

            g_path_comp = [];
            % path complementarities
            % TODO: do this cleaner
            if (model.g_comp_path_constraint &&...
                (obj.fe_idx == dims.N_finite_elements(obj.ctrl_idx) || settings.g_path_at_fe))
                g_path_comp = vertcat(g_path_comp, model.g_comp_path_fun(obj.prev_fe.x{end}, obj.u, p_stage, model.v_global));
            end
            for j=1:dims.n_s-settings.right_boundary_point_explicit
                % TODO: there has to be a better way to do this.
                if model.g_comp_path_constraint && settings.g_path_at_stg
                    g_path_comp = vertcat(g_path_comp, model.g_comp_path_fun(obj.x{j}, obj.u, p_stage, model.v_global));
                end
            end

            g_cross_comp = [];
            % TODO Implement other modes
            if ~settings.use_fesd
                % TODO vectorize?
                for j=1:dims.n_s
                    theta = vertcat(obj.theta{j,:});
                    lambda = vertcat(obj.lambda{j,:});
                    g_cross_comp = vertcat(g_cross_comp, diag(theta)*lambda);
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
                for r=1:dims.n_sys
                    sum_theta = obj.sumTheta(r);
                    sum_lambda = obj.sumLambda(r)
                    g_cross_comp = vertcat(g_cross_comp, diag(sum_theta)*sum_lambda);
                end

            elseif settings.cross_comp_mode == 8
                for r=1:dims.n_sys
                    sum_theta = obj.sumTheta(r);
                    sum_lambda = obj.sumLambda(r);
                    g_cross_comp = vertcat(g_cross_comp, dot(sum_theta, sum_lambda));
                end
            elseif settings.cross_comp_mode > 8
                % do nothing in this case
                return
            end

            g_comp = vertcat(g_cross_comp, g_path_comp);

            n_cross_comp = length(g_cross_comp);
            n_path_comp = length(g_path_comp);
            n_comp = n_cross_comp + n_path_comp;
            %
            if ismember(settings.mpcc_mode, MpccMode.elastic_ell_1)
                s_elastic = define_casadi_symbolic(settings.casadi_symbolic_mode,...
                    ['s_elastic_' num2str(obj.ctrl_idx) '_' num2str(obj.fe_idx)],...
                    n_comp);
                obj.addVariable(s_elastic,...
                    'elastic',...
                    settings.s_elastic_min*ones(n_comp,1),...
                    settings.s_elastic_max*ones(n_comp,1),...
                    settings.s_elastic_0*ones(n_comp,1));
            end

            [g_comp, g_comp_lb, g_comp_ub, cost] = reformulate_complementarities(g_comp, settings.mpcc_mode, sigma_p, s_elastic);

            % Add reformulated constraints
            obj.addConstraint(g_comp, g_comp_lb, g_comp_ub);

            % If We need to add a cost from the reformulation do that as needed;
            if settings.mpcc_mode == MpccMode.ell_1_penalty
                if settings.objective_scaling_direct
                    obj.cost = obj.cost + (1/sigma_p)*cost;
                else
                    obj.cost = sigma_p*obj.cost + cost;
                end
            end
        end

        function stepEquilibration(obj, sigma_p, rho_h_p)
            import casadi.*
            model = obj.model;
            settings = obj.settings;
            dims = obj.dims;

            % TODO implement other modes!
            if settings.use_fesd && obj.fe_idx > 1
                nu = obj.nu_vector;
                delta_h_ki = obj.h - obj.prev_fe.h;
                if settings.step_equilibration == StepEquilibrationMode.heuristic_mean
                    h_fe = model.T / (sum(dims.N_finite_elements)); % TODO this may be a bad idea if using different N_fe. may want to issue warning in that case
                    obj.cost = obj.cost + rho_h_p * (obj.h - h_fe).^2;
                elseif settings.step_equilibration ==  StepEquilibrationMode.heuristic_diff
                    obj.cost = obj.cost + rho_h_p * delta_h_ki.^2;
                elseif settings.step_equilibration == StepEquilibrationMode.l2_relaxed_scaled
                    obj.cost = obj.cost + rho_h_p * tanh(nu/settings.step_equilibration_sigma) * delta_h_ki.^2;
                elseif settings.step_equilibration == StepEquilibrationMode.l2_relaxed
                    obj.cost = obj.cost + rho_h_p * nu * delta_h_ki.^2
                elseif settings.step_equilibration == StepEquilibrationMode.direct
                    obj.addConstraint(nu*delta_h_ki, 0, 0);
                elseif settings.step_equilibration == StepEquilibrationMode.direct_homotopy
                    obj.addConstraint([nu*delta_h_ki-sigma_p;-nu*delta_h_ki-sigma_p],...
                        [-inf;-inf],...
                        [0;0]);
                elseif settings.step_equilibration == StepEquilibrationMode.direct_homotopy_lift
                    obj.addConstraint([obj.nu_lift-nu;obj.nu_lift*delta_h_ki-sigma_p;-obj.nu_lift*delta_h_ki-sigma_p],...
                        [0;-inf;-inf],...
                        [0;0;0]);
                else
                    error("Step equilibration mode not implemented");
                end
            end
        end
    end
end

