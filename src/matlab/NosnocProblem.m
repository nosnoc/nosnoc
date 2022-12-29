classdef NosnocProblem < NosnocFormulationObject
    properties
        % Index vectors
        ind_x
        ind_u
        ind_v
        ind_theta
        ind_lam
        ind_mu
        ind_alpha
        ind_lambda_n
        ind_lambda_p
        ind_nu_lift
        ind_h
        ind_elastic
        ind_sot % index for speed of time variable
        ind_boundary % index of bundary value lambda and mu
        ind_t_final % Time-optimal problems: define auxilairy variable for the final time.

        model
        settings
        dims
        ocp

        sigma_p

        p

        fe0
        stages
    end

    properties(Dependent, SetAccess=private, Hidden)
        u
    end
    methods
        function obj = NosnocProblem(settings, dims, model)
            import casadi.*
            obj@NosnocFormulationObject();

            obj.ind_x = {};
            obj.ind_v = {};
            obj.ind_theta = {};
            obj.ind_lam = {};
            obj.ind_mu = {};
            obj.ind_alpha = {};
            obj.ind_lambda_n = {};
            obj.ind_lambda_p = {};
            obj.ind_h = {};
            obj.ind_nu_lift = {};
            obj.ind_sot = {};

            obj.settings = settings;
            obj.model = model;
            obj.dims = dims;
            obj.ocp = []; % TODO create ocp objects

            obj.stages = {};

            sigma_p = SX.sym('sigma_p');
            obj.p = sigma_p;

            obj.createPrimalVariables();

            for ii=1:dims.N_stages
                stage = obj.stages{ii};
                Uk = obj.u{ii};
                for fe=stage
                    % TODO: OCP
                    % 1) Stewart Runge-Kutta discretization
                    fe.forward_simulation(obj.ocp, Uk)

                    % 2) Complementarity Constraints
                    fe.create_complementarity_constraints(sigma_p)

                    % 3) Step Equilibration
                    fe.step_equilibration()

                    % 4) add cost and constraints from FE to problem
                    obj.cost = obj.cost + fe.cost
                    self.add_constraint(fe.g, fe.lbg, fe.ubg)
                end
                if settings.use_fesd && settings.equidistant_control_grid
                    self.addConstraint(sum(arrayfun(@(fe) fe.h, stage)) - h_ctrl_stage)
                end
            end
        end

        % TODO this should be private
        function createPrimalVariables(obj)
            fe0 = FiniteElementZero(obj.settings, obj.model);
            obj.fe0 = fe0;
            
            obj.p = vertcat(obj.p, fe0.lambda);

            self.addVariable(fe0.x{1},...
                             'x',...
                             fe0.lbw(fe0.ind_x{1}),...
                             fe0.ubw(fe0.ind_x{1}),...
                             fe0.w0(fe0.ind_x{1}));

            prev_fe = fe0;
            for ii=1:obj.settings.N_stages
                stage = obj.createControlStage(obj);
                obj.stages = [obj.stages{:}, stage];
                prev_fe = stage{end};
            end
        end

        % TODO this should be private
        function createControlStage(obj, ctrl_idx, prev_fe)
            Uk = SX.sym(['U_' num2str(ctrl_idx)], obj.dims.n_u);
            obj.addVariable(Uk, 'u', obj.model.lbu, obj.model.ubu,...
                            zeros(obj.dims.n_u,1));

            control_stage = {};
            for ii=1:obj.dims.N_finite_elements
                fe = FiniteElement(prev_fe, obj.settings, obj.model, obj.dims, ctrl_idx, ii);
                obj.addFiniteElement(fe);
                control_stage = [control_stage{:}, fe];
                prev_fe = fe;
            end
        end

        % TODO this should be private
        function addFiniteElement(obj, fe)
            w_len = length(obj.w);

            obj.addPrimalVector(fe.w, fe.lbw, fe.ubw, fe.w0);

            obj.ind_h = [obj.ind_h, fe.ind_h+w_len];
            obj.ind_x = [obj.ind_x, increment_indices(fe.ind_x, w_len)];
            obj.ind_v = [obj.ind_v, increment_indices(fe.ind_v, w_len)];
            obj.ind_theta = [obj.ind_theta, increment_indices(fe.ind_theta, w_len)];
            obj.ind_lam = [obj.ind_lam, increment_indices(fe.ind_lam, w_len)];
            obj.ind_mu = [obj.ind_mu, increment_indices(fe.ind_mu, w_len)];
            obj.ind_alpha = [obj.ind_alpha, increment_indices(fe.ind_alpha, w_len)];
            obj.ind_lambda_n = [obj.ind_lambda_n, increment_indices(fe.ind_lambdla_n, w_len)];
            obj.ind_lambda_p = [obj.ind_lambda_p, increment_indices(fe.ind_lambdla_p, w_len)];
            obj.ind_nu_lift = [obj.ind_x, fe.ind_nu_lift+w_len];
        end

        function addPrimalVector(obj, symbolic, lb, ub, initial)
            obj.w = vertcat(obj.w, symbolic);
            obj.lbw = vertcat(obj.lbw, lb);
            obj.lbw = vertcat(obj.ubw, ub);
            obj.lbw = vertcat(obj.w0, initial);
        end
        
        function u = get.u(obj)
            u = cellfun(@(x) obj.w(u), obj.ind_u, 'UniformOutput', false);
        end
    end
end

