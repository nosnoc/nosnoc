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
        sot
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

            if ismember(settings.mpcc_mode, MpccMode.elastic)
                s_elastic = SX.sym('s_elastic',1);
            else
                s_elastic = [];
            end
            obj.createPrimalVariables();

            for ii=1:dims.N_stages
                stage = obj.stages{ii};
                Uk = obj.u{ii};
                if settings.time_rescaling && settings.use_speed_of_time_variables
                    if settings.local_speed_of_time_variable
                        s_sot = obj.sot(ii);
                    else
                        s_sot = obj.sot(1);
                    end
                else
                    s_sot = 1;
                end
                for fe=stage
                    % TODO: OCP
                    % TODO: add path constraints
                    % 1) Stewart Runge-Kutta discretization
                    fe.forwardSimulation(obj.ocp, Uk, s_sot);

                    % 2) Complementarity Constraints
                    fe.createComplementarityConstraints(sigma_p, s_elastic);

                    % 3) Step Equilibration
                    fe.stepEquilibration();

                    % 4) add cost and constraints from FE to problem
                    obj.cost = obj.cost + fe.cost;
                    obj.addConstraint(fe.g, fe.lbg, fe.ubg);
                end
                if settings.use_fesd && settings.equidistant_control_grid
                    obj.addConstraint(sum(vertcat(stage.h)) - model.T / dims.N_stages);
                end
            end

            % Process terminal costs
            try
                last_fe = obj.stages{end}(end);
                obj.cost = obj.cost + model.f_q_T_fun(last_fe.x(end));
            catch
                warning('Terminal cost not defined');
            end
            
            % Process elastic costs
            if ismember(settings.mpcc_mode, MpccMode.elastic)
                obj.addVariable(s_elastic, 'elastic', settings.s_elastic_0, settings.s_elastic_min, settings.s_elastic_max);
                if settings.objective_scaling_direct
                    obj.cost = obj.cost + (1/sigma_p)*s_elastic;
                else
                    obj.cost = sigma_p*obj.cost + s_elastic;
                end
            end
            if ismember(settings.mpcc_mode, MpccMode.elastic)
                sum_s_elastic = [] % TODO
                obj.addVariable(s_elastic, 'elastic', settings.s_elastic_0, settings.s_elastic_min, settings.s_elastic_max);
                if settings.objective_scaling_direct
                    obj.cost = obj.cost + (1/sigma_p)*sum_s_elastic;
                else
                    obj.cost = sigma_p*obj.cost + sum_s_elastic;
                end
            end
        end

        % TODO this should be private
        function createPrimalVariables(obj)
            fe0 = FiniteElementZero(obj.settings, obj.dims, obj.model);
            obj.fe0 = fe0;
            
            obj.p = vertcat(obj.p, fe0.lambda{1,:});

            obj.addVariable(fe0.x{1},...
                            'x',...
                            fe0.lbw(fe0.ind_x{1}),...
                            fe0.ubw(fe0.ind_x{1}),...
                            fe0.w0(fe0.ind_x{1}));

            prev_fe = fe0;
            for ii=1:obj.dims.N_stages
                stage = obj.createControlStage(ii, prev_fe);
                obj.stages = [obj.stages{:}, {stage}];
                prev_fe = stage(end);
            end
        end

        % TODO this should be private
        function control_stage = createControlStage(obj, ctrl_idx, prev_fe)
            import casadi.*
            Uk = SX.sym(['U_' num2str(ctrl_idx)], obj.dims.n_u);
            obj.addVariable(Uk, 'u', obj.model.lbu, obj.model.ubu,...
                            zeros(obj.dims.n_u,1));

            if obj.settings.time_rescaling && obj.settings.use_speed_of_time_variables
                if obj.settings.local_speed_of_time_variable
                    % at every stage
                    s_sot_k = SX.sym(['s_sot_' num2str(ctrl_idx)], 1);
                    obj.addVariable(s_sot_k,...
                                    'sot',...
                                    obj.settings.s_sot_min,...
                                    obj.settings.s_sot_max,...
                                    obj.settings.s_sot0,...
                                    ctrl_idx);
                    obj.cost = obj.cost + (s_sot_k-1)^2;
                else
                    if ctrl_idx == 0
                        % only once
                        s_sot = SX.sym('s_sot', 1);
                        obj.addVariable(s_sot,...
                                        'sot',...
                                        obj.settings.s_sot_min,...
                                        obj.settings.s_sot_max,...
                                        obj.settings.s_sot0);
                        obj.cost = obj.cost + (s_sot-1)^2;
                    end
                end
            end
            
            control_stage = [];
            for ii=1:obj.dims.N_finite_elements
                fe = FiniteElement(prev_fe, obj.settings, obj.model, obj.dims, ctrl_idx, ii);
                obj.addFiniteElement(fe);
                control_stage = [control_stage, fe];
                prev_fe = fe;
            end
        end

        % TODO this should be private
        function addFiniteElement(obj, fe)
            w_len = length(obj.w);

            obj.addPrimalVector(fe.w, fe.lbw, fe.ubw, fe.w0);

            obj.ind_h = [obj.ind_h, fe.ind_h+w_len];
            obj.ind_x = [obj.ind_x; increment_indices(fe.ind_x, w_len)];
            obj.ind_v = [obj.ind_v, increment_indices(fe.ind_v, w_len)];
            obj.ind_theta = [obj.ind_theta, increment_indices(fe.ind_theta, w_len)];
            obj.ind_lam = [obj.ind_lam, increment_indices(fe.ind_lam, w_len)];
            obj.ind_mu = [obj.ind_mu, increment_indices(fe.ind_mu, w_len)];
            obj.ind_alpha = [obj.ind_alpha, increment_indices(fe.ind_alpha, w_len)];
            obj.ind_lambda_n = [obj.ind_lambda_n, increment_indices(fe.ind_lambda_n, w_len)];
            obj.ind_lambda_p = [obj.ind_lambda_p, increment_indices(fe.ind_lambda_p, w_len)];
            obj.ind_nu_lift = [obj.ind_x, fe.ind_nu_lift+w_len];
        end

        function addPrimalVector(obj, symbolic, lb, ub, initial)
            obj.w = vertcat(obj.w, symbolic);
            obj.lbw = vertcat(obj.lbw, lb);
            obj.ubw = vertcat(obj.ubw, ub);
            obj.w0 = vertcat(obj.w0, initial);
        end
        
        function u = get.u(obj)
            u = obj.w(obj.ind_u);
        end

        function sot = get.sot(obj)
            sot = obj.w(obj.ind_sot);
        end
        
        function print(obj)
            disp("g");
            print_casadi_vector(obj.g);
            disp('lbg, ubg');
            disp([length(obj.lbg), length(obj.ubg)]);
            disp([obj.lbg, obj.ubg]);

            disp("w");
            print_casadi_vector(obj.w);
            disp('lbw, ubw');
            disp([obj.lbw, obj.ubw]);
            disp('w0');
            disp(obj.w0);

            disp('objective');
            disp(obj.cost);
        end
    end
end

