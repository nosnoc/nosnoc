classdef NosnocProblem < NosnocFormulationObject
    properties
        % Index vectors
        ind_x
        ind_u
        ind_v
        ind_z
        ind_h
        ind_elastic
        ind_sot % index for speed of time variable
        ind_boundary % index of bundary value lambda and mu
        ind_t_final % Time-optimal problems: define auxilairy variable for the final time.

        model
        settings
        dims

        sigma_p

        p
        
        stages
    end

    properties(Dependent, SetAccess=private, Hidden)
        u
    end
    methods
        function obj = NosnocProblem()
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
            
        end

        % TODO this should be private
        function createControlStage(obj)
            
        end

        % TODO this should be private
        function addFiniteElement(obj, fe)
            
        end

        function u = get.u(obj)
            u = cellfun(@(x) obj.w(u), obj.ind_u, 'UniformOutput', false);
        end
    end
end

