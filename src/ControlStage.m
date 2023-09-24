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
        ind_L_vn
        ind_P_vt
        ind_N_vt
        ind_Alpha_vt
        % misc
        ind_h
        ind_nu_lift
        ind_elastic
        ind_sot
        ind_comp_lift

        % Index of this control stage.
        ctrl_idx

        % Control stage g.
        ind_stage
        
        % Problem data
        model
        problem_options
        dims
        ocp
    end

    methods
        % TODO: This probably should take less arguments somehow. Maybe a store of "global_variables" to be
        % added along with v_global. This will likely be done when I get around to cleaning up the `model` struct.
        function obj = ControlStage(prev_fe, problem_options, model, dims, ctrl_idx, s_sot, T_final, rho_h_p, rho_sot_p)
            import casadi.*
            obj@NosnocFormulationObject();
            
            obj.problem_options = problem_options;
            obj.model = model;
            obj.dims = dims;

            obj.ctrl_idx = ctrl_idx;
            
            obj.Uk = define_casadi_symbolic(problem_options.casadi_symbolic_mode, ['U_' num2str(ctrl_idx)], obj.dims.n_u);
            obj.addVariable(obj.Uk, 'u', obj.model.lbu, obj.model.ubu,...
                            zeros(obj.dims.n_u,1));

            p_stage = vertcat(model.p_global, model.p_time_var_stages(:, ctrl_idx));
            if problem_options.time_rescaling && problem_options.use_speed_of_time_variables
                if problem_options.local_speed_of_time_variable
                    % at every stage
                    s_sot = define_casadi_symbolic(problem_options.casadi_symbolic_mode, ['s_sot_' num2str(ctrl_idx)], 1);
                    obj.addVariable(s_sot,...
                                    'sot',...
                                    problem_options.s_sot_min,...
                                    problem_options.s_sot_max,...
                                    problem_options.s_sot0);
                    if problem_options.time_freezing
                        obj.augmented_objective = obj.augmented_objective + rho_sot_p*(s_sot-1)^2;
                    end
                end
            else
                s_sot = 1;
            end
            
            obj.stage = [];
            for ii=1:obj.problem_options.N_finite_elements
                fe = FiniteElement(prev_fe, obj.problem_options, obj.model, obj.dims, ctrl_idx, ii, T_final);
                % 1) Runge-Kutta discretization
                fe.forwardSimulation(obj.ocp, obj.Uk, s_sot, p_stage);

                % 2) Complementarity Constraints
                fe.createComplementarityConstraints(p_stage);

                % 3) Step Equilibration
                fe.stepEquilibration(rho_h_p);

                % 4) add finite element variables
                obj.addFiniteElement(fe);
                
                % 5) add augmented_objective, objective and, constraints from FE to problem
                obj.augmented_objective = obj.augmented_objective + fe.augmented_objective;
                obj.objective = obj.objective + fe.objective;
                
                obj.addConstraint(fe.g, fe.lbg, fe.ubg);
                
                obj.stage = [obj.stage, fe];
                prev_fe = fe;
            end

            % least squares cost
            obj.augmented_objective = obj.augmented_objective + (problem_options.T/problem_options.N_stages)*model.f_lsq_x_fun(obj.stage(end).x{end},model.x_ref_val(:,obj.ctrl_idx), p_stage);
            obj.objective = obj.objective + (problem_options.T/problem_options.N_stages)*model.f_lsq_x_fun(obj.stage(end).x{end},model.x_ref_val(:,obj.ctrl_idx), p_stage);
            if dims.n_u > 0
                obj.augmented_objective = obj.augmented_objective + (problem_options.T/problem_options.N_stages)*model.f_lsq_u_fun(obj.Uk,model.u_ref_val(:,obj.ctrl_idx), p_stage);
                obj.objective = obj.objective + (problem_options.T/problem_options.N_stages)*model.f_lsq_u_fun(obj.Uk,model.u_ref_val(:,obj.ctrl_idx), p_stage);
            end
            
            % TODO: combine this into a function
            if problem_options.use_fesd && problem_options.equidistant_control_grid
                if ~problem_options.time_optimal_problem
                    obj.addConstraint(sum(vertcat(obj.stage.h)) - problem_options.h, 'type', 'stage');
                elseif ~problem_options.time_freezing
                    if problem_options.use_speed_of_time_variables
                        obj.addConstraint(sum(vertcat(obj.stage.h)) - problem_options.h, 'type', 'stage')
                        obj.addConstraint(sum(s_sot*vertcat(obj.stage.h)) - T_final/problem_options.N_stages, 'type', 'stage');
                    else
                        obj.addConstraint(sum(vertcat(obj.stage.h)) - T_final/problem_options.N_stages, 'type', 'stage');
                    end
                end
            end
            if problem_options.time_freezing && problem_options.stagewise_clock_constraint
                if problem_options.time_optimal_problem
                    obj.addConstraint(fe.x{end}(end) - ctrl_idx*(T_final/problem_options.N_stages) + model.x0(end), 'type', 'stage');
                else
                    obj.addConstraint(fe.x{end}(end) - ctrl_idx*problem_options.h + model.x0(end), 'type', 'stage');
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
            obj.ind_L_vn = [obj.ind_L_vn; transpose(flatten_sys(increment_indices(fe.ind_L_vn, w_len)))];
            obj.ind_P_vt = [obj.ind_P_vt; transpose(flatten_sys(increment_indices(fe.ind_P_vt, w_len)))];
            obj.ind_N_vt = [obj.ind_N_vt; transpose(flatten_sys(increment_indices(fe.ind_N_vt, w_len)))];
            obj.ind_Alpha_vt = [obj.ind_Alpha_vt; transpose(flatten_sys(increment_indices(fe.ind_Alpha_vt, w_len)))];
            
            obj.ind_nu_lift = [obj.ind_nu_lift, {fe.ind_nu_lift+w_len}];
            obj.ind_comp_lift = [obj.ind_nu_lift, {fe.ind_comp_lift+w_len}];
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
        
        function [u, lbu, ubu, u0] = u(obj)
            u = obj.w(obj.ind_u);
            lbu = obj.lbw(obj.ind_u);
            ubu = obj.ubw(obj.ind_u);
            u0 = obj.w0(obj.ind_u);
        end

        function [sot, lbsot, ubsot, sot0] = sot(obj)
            sot = obj.w(obj.ind_sot);
            lbsot = obj.lbw(obj.ind_sot);
            ubsot = obj.ubw(obj.ind_sot);
            sot0 = obj.w0(obj.ind_sot);
        end

        function [g_stage, lbg_stage, ubg_stage] = g_stage(obj)
            g_stage = obj.g(obj.ind_stage);
            lbg_stage = obj.lbg(obj.ind_stage);
            ubg_stage = obj.ubg(obj.ind_stage);
        end

        function json = jsonencode(obj, varargin)
            import casadi.*
            stage_struct = struct(obj);

            stage_struct = rmfield(stage_struct, 'model');
            stage_struct = rmfield(stage_struct, 'dims');
            stage_struct = rmfield(stage_struct, 'problem_options');
            json = jsonencode(stage_struct);
        end
    end
end

