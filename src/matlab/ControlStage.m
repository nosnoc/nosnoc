classdef ControlStage < NosnocFormulationObject 
    properties
        stage
        Uk

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
        ind_gamma
        ind_u
        ind_sot

        ctrl_idx
        
        model
        settings
        dims
        ocp
    end

    methods
        function obj = ControlStage(prev_fe, settings, model, dims, ctrl_idx, s_sot, T_final, sigma_p, rho_h_p, rho_sot_p, s_elastic)
            import casadi.*
            obj@NosnocFormulationObject();
            
            obj.settings = settings;
            obj.model = model;
            obj.dims = dims;

            obj.ctrl_idx = ctrl_idx;
            
            obj.Uk = SX.sym(['U_' num2str(ctrl_idx)], obj.dims.n_u);
            obj.addVariable(obj.Uk, 'u', obj.model.lbu, obj.model.ubu,...
                            zeros(obj.dims.n_u,1));

            if settings.time_rescaling && settings.use_speed_of_time_variables
                if settings.local_speed_of_time_variable
                    % at every stage
                    s_sot = SX.sym(['s_sot_' num2str(ctrl_idx)], 1);
                    obj.addVariable(s_sot,...
                                    'sot',...
                                    settings.s_sot_min,...
                                    settings.s_sot_max,...
                                    settings.s_sot0);
                    if settings.time_freezing
                        obj.cost = obj.cost + obj.rho_sot_p*(s_sot_k-1)^2;
                    end
                end
            else
                s_sot = 1;
            end
            
            obj.stage = [];
            for ii=1:obj.dims.N_finite_elements
                fe = FiniteElement(prev_fe, obj.settings, obj.model, obj.dims, ctrl_idx, ii, T_final);
                % TODO: OCP
                % 1) Stewart Runge-Kutta discretization
                fe.forwardSimulation(obj.ocp, obj.Uk, s_sot);

                % 2) Complementarity Constraints
                fe.createComplementarityConstraints(sigma_p, s_elastic);

                % 3) Step Equilibration
                fe.stepEquilibration(sigma_p, rho_h_p);
                % 4) add finite element variables
                obj.addFiniteElement(fe);
                
                % 5) add cost and constraints from FE to problem
                obj.cost = obj.cost + fe.cost;
                obj.addConstraint(fe.g, fe.lbg, fe.ubg);
                
                obj.stage = [obj.stage, fe];
                prev_fe = fe;
            end

            obj.createComplementarityConstraints(sigma_p, s_elastic);

            % least squares cost
            obj.cost = obj.cost + (model.T/dims.N_stages)*model.f_lsq_x_fun(obj.stage(end).x{end},model.x_ref_val(:,obj.ctrl_idx));
            if dims.n_u > 0
                obj.cost = obj.cost + (model.T/dims.N_stages)*model.f_lsq_u_fun(obj.Uk,model.u_ref_val(:,ii));
            end
            
            % TODO: combine this into a function
            if settings.use_fesd && settings.equidistant_control_grid
                if ~settings.time_optimal_problem
                    obj.addConstraint(sum(vertcat(obj.stage.h)) - model.h);
                elseif ~settings.time_freezing
                    if settings.use_speed_of_time_variables
                        obj.addConstraint(sum(vertcat(obj.stage.h)) - model.h)
                        obj.addConstraint(sum(s_sot*vertcat(obj.stage.h)) - T_final/dims.N_stages);
                    else
                        obj.addConstraint(sum(vertcat(obj.stage.h)) - T_final/dims.N_stages);
                    end
                end
            end
            if settings.time_freezing && settings.stagewise_clock_constraint
                if time_optimal_problem
                    obj.addConstraint(fe.x{end}(end) - ctrl_idx*(T_final/dims.N_stages) + model.x0(end));
                else
                    obj.addConstraint(fe.x{end}(end) - ctrl_idx*model.h + model.x0(end));
                end
            end
        end 

        % TODO this should be private
        function addFiniteElement(obj, fe)
            w_len = size(obj.w, 1);

            obj.addPrimalVector(fe.w, fe.lbw, fe.ubw, fe.w0);

            obj.ind_h = [obj.ind_h, fe.ind_h+w_len];
            obj.ind_x = [obj.ind_x, increment_indices(fe.ind_x, w_len)];
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
            lens = [size(symbolic,1), size(lb,1), size(ub,1), size(initial,1)];
            if ~all(lens == lens(1))
                symbolic
                lb
                ub
                initial
                error("mismatched dims")
            end
            obj.w = vertcat(obj.w, symbolic);
            obj.lbw = vertcat(obj.lbw, lb);
            obj.ubw = vertcat(obj.ubw, ub);
            obj.w0 = vertcat(obj.w0, initial);
        end
        
        function createComplementarityConstraints(obj, sigma_p, s_elastic)
            import casadi.*           
            model = obj.model;
            settings = obj.settings;
            dims = obj.dims;

            % path constraints handled at fe level not here
            
            g_cross_comp = SX([]);
            % TODO Implement other modes
            if ~settings.use_fesd || settings.cross_comp_mode < 9 || settings.cross_comp_mode > 10
                % Do nothing, handled at the FE level
                % along with modes 1-8
                return
                
            elseif settings.cross_comp_mode == 9
                for r=1:dims.n_sys
                    g_r = 0;
                    for fe=obj.stage
                        g_r = g_r + diag(fe.sumTheta(r))*fe.sumLambda(r);
                    end
                    g_cross_comp = vertcat(g_cross_comp, g_r);
                end
            elseif settings.cross_comp_mode == 10
                for r=1:dims.n_sys
                    g_r = 0;
                    for fe=obj.stage
                        g_r = g_r + dot(fe.sumTheta(r),fe.sumLambda(r));
                    end
                    g_cross_comp = vertcat(g_cross_comp, g_r);
                end
            elseif settings.cross_comp_mode == 11
                error('TODO: not implemented');
            elseif settings.cross_comp_mode == 12
                error('TODO: not implemented');
            end

            g_comp = g_cross_comp;
            n_comp = length(g_cross_comp);
            
            %
            if ismember(settings.mpcc_mode, MpccMode.elastic_ell_1)
                s_elastic = SX.sym(['s_elastic_' num2str(obj.ctrl_idx) '_' num2str(obj.fe_idx)], n_comp);
                obj.addVariable(s_elastic,...
                                'elastic',...
                                settings.s_elastic_min*ones(n_comp,1),...
                                settings.s_elastic_max*ones(n_comp,1),...
                                settings.s_elastic_0*ones(n_comp,1));
            end
            
            % Do MPCC formulation
            % TODO this should be done on the problem level, and take into account passed in vars. fine for now.
            if settings.mpcc_mode == MpccMode.Scholtes_ineq
                g_comp = g_comp - sigma_p;
                g_comp_ub = zeros(n_comp,1);
                g_comp_lb = -inf * ones(n_comp,1);
            elseif settings.mpcc_mode == MpccMode.Scholtes_eq
                g_comp = g_comp - sigma_p;
                g_comp_ub = zeros(n_comp,1);
                g_comp_lb = zeros(n_comp,1);
            elseif settings.mpcc_mode == MpccMode.ell_1_penalty
                if settings.objective_scaling_direct
                    obj.cost = obj.cost + (1/sigma_p)*sum(g_comp);
                else
                    obj.cost = sigma_p*obj.cost + sum(g_comp);
                end
            elseif settings.mpcc_mode == MpccMode.elastic_ineq
                g_comp = g_comp - s_elastic*ones(n_comp,1);
                g_comp_ub = zeros(n_comp,1);
                g_comp_lb = -inf * ones(n_comp,1);
            elseif settings.mpcc_mode == MpccMode.elastic_eq
                g_comp = g_comp - s_elastic*ones(n_comp,1);
                g_comp_ub = zeros(n_comp,1);
                g_comp_lb = zeros(n_comp,1);
            elseif settings.mpcc_mode == MpccMode.elastic_two_sided
                g_comp = [g_comp-s_elastic*ones(n_comp,1);g_comp+s_elastic*ones(n_comp,1)];
                g_comp_ub = [zeros(n_comp,1); inf*ones(n_comp,1)];
                g_comp_lb = [-inf*ones(n_comp,1);  zeros(n_comp,1)];
                % TODO these should be merged probably
            elseif settings.mpcc_mode == MpccMode.elastic_ell_1_ineq
                g_comp = g_comp - s_elastic;
                g_comp_ub = zeros(n_comp,1);
                g_comp_lb = -inf * ones(n_comp,1);
            elseif settings.mpcc_mode == MpccMode.elastic_ell_1_eq
                g_comp = g_comp - s_elastic;
                g_comp_ub = zeros(n_comp,1);
                g_comp_lb = zeros(n_comp,1);
            elseif settings.mpcc_mode == MpccMode.elastic_ell_1_two_sided
                g_comp = [g_comp-s_elastic;g_comp+s_elastic];
                g_comp_ub = [zeros(n_comp,1); inf*ones(n_comp,1)];
                g_comp_lb = [-inf*ones(n_comp,1);  zeros(n_comp,1)];
            end

            if settings.mpcc_mode ~= MpccMode.ell_1_penalty
                obj.addConstraint(g_comp, g_comp_lb, g_comp_ub);
            end
        end
    end
end

