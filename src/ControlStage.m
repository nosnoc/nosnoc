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
classdef ControlStage < NosnocFormulationObject
    properties
        stage
        Uk
        % Index vectors
        ind_x
        ind_v
        ind_z
        ind_u
        % Specific to Stewarts representation
        ind_theta
        ind_lam
        ind_mu
        % Specific to Step representation
        ind_alpha
        ind_lambda_n
        ind_lambda_p
        ind_beta
        ind_theta_step
        % Speficif to CLS representation
        ind_lambda_normal
        ind_lambda_tangent
        ind_y_gap
        % friction multipliers and lifting
        % conic
        ind_gamma
        ind_beta_conic
        % poly
        ind_gamma_d
        ind_beta_d
        ind_delta_d
        % variables related to conic
        ind_p_vt
        ind_n_vt
        ind_alpha_vt
        % variables only at element boundary
        ind_x_left_bp
        ind_Y_gap
        ind_Lambda_normal
        ind_Lambda_tangent
        ind_Gamma
        ind_Gamma_d
        ind_Beta_conic
        ind_Beta_d
        ind_Delta_d
%         ind_P_vn
%         ind_N_vn
        ind_L_vn
        ind_P_vt
        ind_N_vt
        ind_Alpha_vt
        % misc
        ind_h
        ind_nu_lift
        ind_elastic
        ind_sot

        % Index of this control stage.
        ctrl_idx

        % Problem data
        model
        settings
        dims
        ocp
    end

    methods
        % TODO: This probably should take less arguments somehow. Maybe a store of "global_variables" to be
        % added along with v_global. This will likely be done when I get around to cleaning up the `model` struct.
        function obj = ControlStage(prev_fe, settings, model, dims, ctrl_idx, s_sot, T_final, sigma_p, rho_h_p, rho_sot_p, s_elastic)
            import casadi.*
            obj@NosnocFormulationObject();
            
            obj.settings = settings;
            obj.model = model;
            obj.dims = dims;

            obj.ctrl_idx = ctrl_idx;
            
            obj.Uk = define_casadi_symbolic(settings.casadi_symbolic_mode, ['U_' num2str(ctrl_idx)], obj.dims.n_u);
            obj.addVariable(obj.Uk, 'u', obj.model.lbu, obj.model.ubu,...
                            zeros(obj.dims.n_u,1));

            p_stage = vertcat(model.p_global, model.p_time_var_stages(:, ctrl_idx));
            if settings.time_rescaling && settings.use_speed_of_time_variables
                if settings.local_speed_of_time_variable
                    % at every stage
                    s_sot = define_casadi_symbolic(settings.casadi_symbolic_mode, ['s_sot_' num2str(ctrl_idx)], 1);
                    obj.addVariable(s_sot,...
                                    'sot',...
                                    settings.s_sot_min,...
                                    settings.s_sot_max,...
                                    settings.s_sot0);
                    if settings.time_freezing
                        obj.cost = obj.cost + rho_sot_p*(s_sot-1)^2;
                    end
                end
            else
                s_sot = 1;
            end
            
            obj.stage = [];
            for ii=1:obj.settings.N_finite_elements
                fe = FiniteElement(prev_fe, obj.settings, obj.model, obj.dims, ctrl_idx, ii, T_final);
                % 1) Runge-Kutta discretization
                fe.forwardSimulation(obj.ocp, obj.Uk, s_sot, p_stage);

                % 2) Complementarity Constraints
                fe.createComplementarityConstraints(sigma_p, s_elastic, p_stage);

                % 3) Step Equilibration
                fe.stepEquilibration(sigma_p, rho_h_p);

                % 4) add finite element variables
                obj.addFiniteElement(fe);
                
                % 5) add cost, objective and, constraints from FE to problem
                obj.cost = obj.cost + fe.cost;
                obj.objective = obj.objective + fe.objective;
                
                obj.addConstraint(fe.g, fe.lbg, fe.ubg);
                
                obj.stage = [obj.stage, fe];
                prev_fe = fe;
            end

            obj.createComplementarityConstraints(sigma_p, s_elastic);

            % least squares cost
            obj.cost = obj.cost + (model.T/settings.N_stages)*model.f_lsq_x_fun(obj.stage(end).x{end},model.x_ref_val(:,obj.ctrl_idx), p_stage);
            obj.objective = obj.objective + (model.T/settings.N_stages)*model.f_lsq_x_fun(obj.stage(end).x{end},model.x_ref_val(:,obj.ctrl_idx), p_stage);
            if dims.n_u > 0
                obj.cost = obj.cost + (model.T/settings.N_stages)*model.f_lsq_u_fun(obj.Uk,model.u_ref_val(:,obj.ctrl_idx), p_stage);
                obj.objective = obj.objective + (model.T/settings.N_stages)*model.f_lsq_u_fun(obj.Uk,model.u_ref_val(:,obj.ctrl_idx), p_stage);
            end
            
            % TODO: combine this into a function
            if settings.use_fesd && settings.equidistant_control_grid
                if ~settings.time_optimal_problem
                    obj.addConstraint(sum(vertcat(obj.stage.h)) - model.h);
                elseif ~settings.time_freezing
                    if settings.use_speed_of_time_variables
                        obj.addConstraint(sum(vertcat(obj.stage.h)) - model.h)
                        obj.addConstraint(sum(s_sot*vertcat(obj.stage.h)) - T_final/settings.N_stages);
                    else
                        obj.addConstraint(sum(vertcat(obj.stage.h)) - T_final/settings.N_stages);
                    end
                end
            end
            if settings.time_freezing && settings.stagewise_clock_constraint
                if settings.time_optimal_problem
                    obj.addConstraint(fe.x{end}(end) - ctrl_idx*(T_final/settings.N_stages) + model.x0(end));
                else
                    obj.addConstraint(fe.x{end}(end) - ctrl_idx*model.h + model.x0(end));
                end
            end
        end

        % TODO this should be private
        function addFiniteElement(obj, fe)
            w_len = size(obj.w, 1);

            obj.addPrimalVector(fe.w, fe.lbw, fe.ubw, fe.w0);

            flatten_sys =  @(ind) arrayfun(@(ii) [ind{ii,:}], transpose(1:size(ind,1)), 'uni', 0);

            obj.ind_h = [obj.ind_h, {fe.ind_h+w_len}];
            obj.ind_x = [obj.ind_x; increment_indices(fe.ind_x.', w_len)];
            obj.ind_v = [obj.ind_v; increment_indices(fe.ind_v.', w_len)];
            obj.ind_theta = [obj.ind_theta; transpose(flatten_sys(increment_indices(fe.ind_theta, w_len)))];
            obj.ind_lam = [obj.ind_lam; transpose(flatten_sys(increment_indices(fe.ind_lam, w_len)))];
            obj.ind_mu = [obj.ind_mu; transpose(flatten_sys(increment_indices(fe.ind_mu, w_len)))];
            obj.ind_alpha = [obj.ind_alpha; transpose(flatten_sys(increment_indices(fe.ind_alpha, w_len)))];
            obj.ind_lambda_n = [obj.ind_lambda_n; transpose(flatten_sys(increment_indices(fe.ind_lambda_n, w_len)))];
            obj.ind_lambda_p = [obj.ind_lambda_p; transpose(flatten_sys(increment_indices(fe.ind_lambda_p, w_len)))];
            obj.ind_beta = [obj.ind_beta; transpose(flatten_sys(increment_indices(fe.ind_beta, w_len)))];
            obj.ind_theta_step = [obj.ind_theta_step; transpose(flatten_sys(increment_indices(fe.ind_theta_step, w_len)))];
            obj.ind_z = [obj.ind_z; transpose(flatten_sys(increment_indices(fe.ind_z, w_len)))];
            obj.ind_lambda_normal = [obj.ind_lambda_normal; transpose(flatten_sys(increment_indices(fe.ind_lambda_normal, w_len)))];
            obj.ind_lambda_tangent = [obj.ind_lambda_tangent; transpose(flatten_sys(increment_indices(fe.ind_lambda_tangent, w_len)))];
            obj.ind_y_gap = [obj.ind_y_gap; transpose(flatten_sys(increment_indices(fe.ind_y_gap, w_len)))];
            obj.ind_gamma = [obj.ind_gamma; transpose(flatten_sys(increment_indices(fe.ind_gamma, w_len)))];
            obj.ind_beta_conic = [obj.ind_beta_conic; transpose(flatten_sys(increment_indices(fe.ind_beta_conic, w_len)))];
            obj.ind_gamma_d = [obj.ind_gamma_d; transpose(flatten_sys(increment_indices(fe.ind_gamma_d, w_len)))];
            obj.ind_beta_d = [obj.ind_beta_d; transpose(flatten_sys(increment_indices(fe.ind_beta_d, w_len)))];
            obj.ind_delta_d = [obj.ind_delta_d; transpose(flatten_sys(increment_indices(fe.ind_delta_d, w_len)))];
            obj.ind_p_vt = [obj.ind_p_vt; transpose(flatten_sys(increment_indices(fe.ind_p_vt, w_len)))];
            obj.ind_n_vt = [obj.ind_n_vt; transpose(flatten_sys(increment_indices(fe.ind_n_vt, w_len)))];
            obj.ind_alpha_vt = [obj.ind_alpha_vt; transpose(flatten_sys(increment_indices(fe.ind_alpha_vt, w_len)))];

            obj.ind_x_left_bp = [obj.ind_x_left_bp; transpose(flatten_sys(increment_indices(fe.ind_x_left_bp, w_len)))];
            obj.ind_Y_gap = [obj.ind_Y_gap; transpose(flatten_sys(increment_indices(fe.ind_Y_gap, w_len)))];
            obj.ind_Lambda_normal = [obj.ind_Lambda_normal; transpose(flatten_sys(increment_indices(fe.ind_Lambda_normal, w_len)))];
            obj.ind_Lambda_tangent = [obj.ind_Lambda_tangent; transpose(flatten_sys(increment_indices(fe.ind_Lambda_tangent, w_len)))];
            obj.ind_Gamma = [obj.ind_Gamma; transpose(flatten_sys(increment_indices(fe.ind_Gamma, w_len)))];
            obj.ind_Beta_conic = [obj.ind_Beta_conic; transpose(flatten_sys(increment_indices(fe.ind_Beta_conic, w_len)))];
            obj.ind_Gamma_d = [obj.ind_Gamma_d; transpose(flatten_sys(increment_indices(fe.ind_Gamma_d, w_len)))];
            obj.ind_Beta_d = [obj.ind_Beta_d; transpose(flatten_sys(increment_indices(fe.ind_Beta_d, w_len)))];
            obj.ind_Delta_d = [obj.ind_Delta_d; transpose(flatten_sys(increment_indices(fe.ind_Delta_d, w_len)))];
%             obj.ind_P_vn = [obj.ind_P_vn; transpose(flatten_sys(increment_indices(fe.ind_P_vn, w_len)))];
%             obj.ind_N_vn = [obj.ind_N_vn; transpose(flatten_sys(increment_indices(fe.ind_N_vn, w_len)))];
            obj.ind_L_vn = [obj.ind_L_vn; transpose(flatten_sys(increment_indices(fe.ind_L_vn, w_len)))];
            obj.ind_P_vt = [obj.ind_P_vt; transpose(flatten_sys(increment_indices(fe.ind_P_vt, w_len)))];
            obj.ind_N_vt = [obj.ind_N_vt; transpose(flatten_sys(increment_indices(fe.ind_N_vt, w_len)))];
            obj.ind_Alpha_vt = [obj.ind_Alpha_vt; transpose(flatten_sys(increment_indices(fe.ind_Alpha_vt, w_len)))];
            
            obj.ind_nu_lift = [obj.ind_nu_lift, {fe.ind_nu_lift+w_len}];
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

            psi_fun = settings.psi_fun;

            % Generate the s_elastic variable if necessary.
            if ismember(settings.mpcc_mode, MpccMode.elastic_ell_1)
                s_elastic = define_casadi_symbolic(settings.casadi_symbolic_mode, ['s_elastic_' num2str(obj.ctrl_idx)], n_comp);
                obj.addVariable(s_elastic,...
                                'elastic',...
                                settings.s_elastic_min*ones(n_comp,1),...
                                settings.s_elastic_max*ones(n_comp,1),...
                                settings.s_elastic_0*ones(n_comp,1));
            end

            if settings.elasticity_mode == ElasticityMode.NONE
                sigma = sigma_p;
            else
                sigma = s_elastic;
            end
            
            g_cross_comp = [];
            if ~settings.use_fesd || settings.cross_comp_mode < 9 || settings.cross_comp_mode > 10
                % Do nothing, handled at the FE level
                % along with modes 1-8
                return                
            elseif settings.cross_comp_mode == 9
                for r=1:dims.n_sys
                    g_r = 0;
                    nz_r = [];
                    for fe=obj.stage
                        pairs = fe.cross_comp_pairs(:, :, r);
                        expr_cell = cellfun(@(pair) apply_psi(pair, psi_fun, sigma), pairs, 'uni', false);
                        nonzeros = cellfun(@(x) vector_is_zero(x), expr_cell, 'uni', 0);
                        if size([expr_cell{:}], 1) == 0
                            exprs= [];
                        elseif settings.relaxation_method == RelaxationMode.TWO_SIDED
                            exprs_p = cellfun(@(c) c(:,1), expr_cell, 'uni', false);
                            exprs_n = cellfun(@(c) c(:,2), expr_cell, 'uni', false);
                            nonzeros_p = cellfun(@(x) vector_is_zero(x), exprs_p, 'uni', 0);
                            nonzeros_n = cellfun(@(x) vector_is_zero(x), exprs_n, 'uni', 0);
                            nonzeros = [sum([nonzeros_p{:}], 2),sum([nonzeros_n{:}], 2)]';
                            exprs = [sum2([exprs_p{:}]),sum2([exprs_n{:}])]';
                            exprs = exprs(:);
                        else
                            nonzeros = sum([nonzeros{:}], 2);
                            exprs = sum2([expr_cell{:}]);
                        end
                        if isempty(nz_r)
                            nz_r = zeros(size(nonzeros));
                        end
                        g_r = g_r + extract_nonzeros_from_vector(exprs);
                        nz_r = nz_r + nonzeros;
                    end
                    g_r = scale_sigma(g_r, sigma, nz_r);
                    g_cross_comp = vertcat(g_cross_comp, g_r);
                end
            elseif settings.cross_comp_mode == 10
                for r=1:dims.n_sys
                    g_r = 0;
                    nz_r = [];
                    for fe=obj.stage
                        pairs = fe.cross_comp_pairs(:, :, r);
                        expr_cell = cellfun(@(pair) apply_psi(pair, psi_fun, sigma), pairs, 'uni', false);
                        nonzeros = cellfun(@(x) vector_is_zero(x), expr_cell, 'uni', 0);
                        if size([expr_cell{:}], 1) == 0
                            exprs= [];
                        elseif settings.relaxation_method == RelaxationMode.TWO_SIDED
                            exprs_p = cellfun(@(c) c(:,1), expr_cell, 'uni', false);
                            exprs_n = cellfun(@(c) c(:,2), expr_cell, 'uni', false);
                            nonzeros_p = cellfun(@(x) vector_is_zero(x), exprs_p, 'uni', 0);
                            nonzeros_n = cellfun(@(x) vector_is_zero(x), exprs_n, 'uni', 0);
                            nonzeros = [sum(sum([nonzeros_p{:}], 2), 1),sum(sum([nonzeros_n{:}], 2), 1)]';
                            exprs = [sum1(sum2([exprs_p{:}])),sum1(sum2([exprs_n{:}]))]';
                            exprs = exprs(:);
                        else
                            nonzeros = sum(sum([nonzeros{:}], 2),1);
                            exprs = sum1(sum2([expr_cell{:}]));
                        end
                        g_r = g_r + extract_nonzeros_from_vector(exprs);
                        if isempty(nz_r)
                            nz_r = zeros(size(nonzeros));
                        end
                        nz_r = nz_r + nonzeros;
                    end
                    g_r = scale_sigma(g_r, sigma, nz_r);
                    g_cross_comp = vertcat(g_cross_comp, g_r);
                end
            end

            g_comp = g_cross_comp;
            n_comp = length(g_cross_comp);

            [g_comp_lb, g_comp_ub, g_comp] = generate_mpcc_relaxation_bounds(g_comp, settings);
            
            obj.addConstraint(g_comp, g_comp_lb, g_comp_ub);

            % TODO handle ell_1_penalty at top level.
            % If We need to add a cost from the reformulation do that.
            if settings.mpcc_mode == MpccMode.ell_1_penalty
                if settings.objective_scaling_direct
                    obj.cost = obj.cost + (1/sigma_p)*cost;
                else
                    obj.cost = sigma_p*obj.cost + cost;
                end
            end
        end
    end
end

