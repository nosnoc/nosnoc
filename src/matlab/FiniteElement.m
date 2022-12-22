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
        ind_elastic
        ind_boundary % index of bundary value lambda and mu, TODO is this even necessary?

        n_stage
        s_sys
        
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
            obj.ind_elastic = [];
            obj.ind_boundary = [];

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
                if settings.irk_representation == IrkRepresentation.DIFFERENTIAL
                    obj.addVariable(SX.sym(['V_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii)'], dims.n_x),...
                                    'v',...
                                    -inf * ones(dims.n_x),...
                                    inf * ones(dims.n_x),...
                                    zeros(dims.n_x),...
                                    ii);
                end
                if settings.irk_representation == IrkRepresentation.INTEGRAL || settings.lift_irk_differential
                    obj.addVariable(SX.sym(['X_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii)], dims.n_x),...
                                    'x',...
                                    -inf * ones(dims.n_x),...
                                    inf * ones(dims.n_x),...
                                    model.x0,...
                                    ii);
                end
                % algebraic variables
                if settings.pss_mode == PssMode.STEWART
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
                                        zeros(dims.n_f_sys(ij)),...
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
                elseif settings.pss_mode == PssMode.STEP
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
                    % add lambda_n
                    for ij = 1:dims.n_sys
                        obj.addVariable(SX.sym(['lambda_n_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)], dims.n_c_sys(ij)),...
                                        'lambda_n',...
                                        zeros(dims.n_c_sys(ij)),...
                                        inf * ones(dims.n_c_sys(ij)),...
                                        settings.initial_lambda * ones(dims.n_c_sys(ij)),...
                                        ii,...
                                        ij);
                    end
                    % add lambda_p
                    for ij = 1:dims.n_sys
                        obj.addVariable(SX.sym(['lambda_p_' num2str(ctrl_idx) '_' num2str(fe_idx) '_' num2str(ii) '_' num2str(ij)],dims.n_c_sys(ij)),...
                                        'lambda_p',...
                                        zeros(dims.n_c_sys(ij)),...
                                        inf * ones(dims.n_c_sys(ij)),...
                                        settings.initial_mu * ones(dims.n_c_sys(ij)),...
                                        ii,...
                                        ij);
                    end
                end
            end
            % Add right boundary points if needed
            if false % TODO: implement create_right_boundary_point
                if settings.pss_mode == PssMode.STEWART
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
                elseif settings.pss_mode == PssMode.STEP
                    % add lambda_n
                    for ij = 1:dims.n_sys
                        obj.addVariable(SX.sym(['lambda_n_' num2str(ctrl_idx) '_' num2str(fe_idx) '_end_' num2str(ij)], dims.n_c_sys(ij)),...
                                        'lambda_n',...
                                        zeros(dims.n_c_sys(ij)),...
                                        inf * ones(dims.n_c_sys(ij)),...
                                        settings.initial_lambda * ones(dims.n_c_sys(ij)),...
                                        dims.n_s,...
                                        ij);
                    end
                    % add lambda_p
                    for ij = 1:dims.n_sys
                        obj.addVariable(SX.sym(['lambda_p_' num2str(ctrl_idx) '_' num2str(fe_idx) '_end_' num2str(ij)], dims.n_c_sys(ij)),...
                                        'lambda_p',...
                                        zeros(dims.n_c_sys(ij)),...
                                        inf * ones(dims.n_c_sys(ij)),...
                                        settings.initial_mu * ones(dims.n_c_sys(ij)),...
                                        dims.n_s,...
                                        ij);
                    end
                end
                % add final X variables
                obj.addVariable(SX.sym(['X_end_' num2str(ctrl_idx) '_' num2str(fe_idx)], dims.n_x),...
                                'x',...
                                -inf * ones(dims.n_x),...
                                inf * ones(dims.n_x),...
                                model.x0,...
                                ii+1);
            end 
        end
        
        function lambda = get.lambda(obj)
            import casadi.*
            grab = @(l, ln, lp) vertcat(obj.w(l), obj.w(ln), obj.w(lp));

            lambda = cellfun(grab, obj.ind_lam, obj.ind_lambda_p, obj.ind_lambda_n, 'UniformOutput', false);
        end

        function theta = get.theta(obj)
            import casadi.*
            grab = @(t, a) vertcat(obj.w(t), obj.w(a), ones(size(a)) - obj.w(a));

            theta = cellfun(grab, obj.ind_theta, obj.ind_alpha, 'UniformOutput', false);
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

        function sum_theta = sumTheta(obj)
            theta = obj.theta;
            thetas = arrayfun(@(row) vertcat(theta{row, :}), 1:size(theta,1), 'UniformOutput', false);

            sum_theta = sum([thetas{:}], 2);
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

        function forwardSimulation(obj)

        end

        function createComplementarityConstraints(obj)

        end

        function stepEquilibration(obj)

        end
    end
end

