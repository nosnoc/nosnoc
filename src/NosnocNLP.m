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


% TODO cleanup steps:
%      - Create primal variables all at once.
%      - Separate sections into separate functions operating on the `problem` struct/class
%      - time variables should probably not just be lumped into the state, for readability.
%      - remove index in symbolic variable defintions and add instructive
%        names, e.g., Uk -> U,  h_ki -> h_fe, X_ki_stages ->  X_rk_stages
%      - provide instructive names for terminal constraint relaxations
%      - provide more instructive names for cross_comp (match python)

classdef NosnocNLP < NosnocFormulationObject
    properties
        ind_elastic

        % Parameter index variables
        ind_p_x0
        ind_p_global
        ind_p_time_var

        % Problem data
        mpcc
        problem_options

        % original initialization
        w0_original

        % Algorithmic parameters
        sigma_p

        % Algorithmic global variables (time independent)
        s_elastic

        % Parameters
        p
        p0

        % Problem cost function
        cost_fun

        % Problem objective function
        objective_fun

        % Problem constraint function
        g_fun

        % switch indicator function
        nu_fun

        nabla_J
        nabla_J_fun
    end
    % remaining list of TODOs
    % TODO: cleanup/add properties (in all components)
    % TODO: Create solver object, which will interact with setting parameters.

    properties(Dependent, SetAccess=private, Hidden)
        % Properties generated on the fly.

        % casadi symbolics/expresions for u, sot, and nu
        u
        sot
        nu_vector
        cc_vector

        % Indices for all algebraic vars in the problem
        ind_z_all
        ind_x_all
    end

    methods
        function obj = NosnocNLP(solver_options, dims, mpcc)
            import casadi.*
            obj@NosnocFormulationObject();

            obj.mpcc = mpcc;
            obj.solver_options = solver_options;
            obj.ocp = []; % TODO create ocp objects
            
            sigma_p = define_casadi_symbolic(settings.casadi_symbolic_mode, 'sigma_p');
            obj.sigma_p = sigma_p;
            obj.p = [sigma_p;mpcc.p];

            obj.cost = mpcc.cost;

            obj.relax_complementarities(sigma_p);

            % Process elastic costs
            if settings.elasticity_mode == ElasticityMode.ELL_INF
                if settings.objective_scaling_direct
                    obj.cost = obj.cost + (1/sigma_p)*obj.s_elastic;
                else
                    obj.cost = sigma_p*obj.cost + obj.s_elastic;
                end
            end
            if settings.elasticity_mode == ElasticityMode.ELL_1
                sum_s_elastic = 0;
                for k=1:settings.N_stages
                    stage=obj.stages(k);
                    for fe=stage.stage
                        sum_s_elastic = sum_s_elastic + fe.sumElastic;
                    end
                end
                if settings.objective_scaling_direct
                    obj.cost = obj.cost + (1/sigma_p)*sum_s_elastic;
                else
                    obj.cost = sigma_p*obj.cost + sum_s_elastic;
                end
            end

            % TODO Figure out if any of these are needed and cleanup how they are calculated
            %obj.comp_res = Function('comp_res', {obj.w, obj.p}, {J_comp});
            obj.comp_std = Function('comp_std', {obj.w, obj.p}, {J_comp_std});
            obj.comp_fesd = Function('comp_fesd', {obj.w, obj.p}, {J_comp_fesd});
            obj.cost_fun = Function('cost_fun', {obj.w, obj.p}, {obj.cost});
            obj.objective_fun = Function('objective_fun', {obj.w, obj.p}, {obj.objective});
            obj.g_fun = Function('g_fun', {obj.w, obj.p}, {obj.g});

            obj.p0 = [settings.sigma_0; settings.rho_sot; settings.rho_h; settings.rho_terminal; model.T];

            if dims.n_p_global > 0
                obj.p0 = [obj.p0; model.p_global_val];
            end

            if dims.n_p_time_var > 0
                obj.p0 = [obj.p0; model.p_time_var_val(:)];
            end
            obj.w0_original = obj.w0;

            % Define CasADi function for the switch indicator function.
            nu_fun = Function('nu_fun', {obj.w,obj.p},{obj.nu_vector});
            obj.nu_fun = nu_fun;
            
            % create CasADi function for objective gradient.
            nabla_J = obj.cost.jacobian(obj.w);
            nabla_J_fun = Function('nabla_J_fun', {obj.w,obj.p},{nabla_J});
            obj.nabla_J = nabla_J;
            obj.nabla_J_fun = nabla_J_fun;
        end

        % TODO this should be private
        function create_primal_variables(obj)
            import casadi.*
            fe0 = FiniteElementZero(obj.settings, obj.dims, obj.model);
            obj.fe0 = fe0;

            %             obj.p = vertcat(obj.p, fe0.x0, fe0.lambda{1,:},fe0.y_gap{1,:},fe0.gamma{1,:},fe0.gamma_d{1,:},fe0.delta_d{1,:},fe0.p_vt{1,:},fe0.n_vt{1,:});
            obj.p = vertcat(obj.p, fe0.x0, fe0.cross_comp_cont_0{1,:}, fe0.cross_comp_cont_1{1,:},fe0.cross_comp_cont_2{1,:});

            X0 = fe0.x{1};
            obj.addVariable(X0,...
                'x0',...
                fe0.lbw(fe0.ind_x{1}),...
                fe0.ubw(fe0.ind_x{1}),...
                fe0.w0(fe0.ind_x{1}));
            obj.addConstraint(fe0.g, fe0.lbg, fe0.ubg);
            prev_fe = fe0;

            s_sot = [];
            if obj.settings.time_rescaling && obj.settings.use_speed_of_time_variables
                if ~obj.settings.local_speed_of_time_variable
                    s_sot = define_casadi_symbolic(obj.settings.casadi_symbolic_mode, 's_sot', 1);
                    obj.addVariable(s_sot,...
                        'sot',...
                        obj.settings.s_sot_min,...
                        obj.settings.s_sot_max,...
                        obj.settings.s_sot0);
                    if obj.settings.time_freezing
                        obj.cost = obj.cost + obj.rho_sot_p*(s_sot-1)^2;
                    end
                end
            end

            if obj.settings.elasticity_mode == ElasticityMode.ELL_INF
                s_elastic = define_casadi_symbolic(obj.settings.casadi_symbolic_mode, 's_elastic',1);
                obj.s_elastic = s_elastic;
                if obj.settings.elastic_scholtes
                    obj.settings.s_elastic_max = inf;
                    obj.addConstraint(s_elastic-obj.sigma_p,-inf,0);
                end
                obj.addVariable(s_elastic, 'elastic', obj.settings.s_elastic_min, obj.settings.s_elastic_max, obj.settings.s_elastic_0);
            else
                s_elastic = [];
            end

            for ii=1:obj.settings.N_stages
                % TODO: maybe this should be a function
                stage = ControlStage(prev_fe, obj.settings, obj.model, obj.dims, ii, s_sot, obj.T_final, obj.sigma_p, obj.rho_h_p, obj.rho_sot_p, s_elastic);
                obj.stages = [obj.stages, stage];

                obj.addControlStage(stage);
                obj.cost = obj.cost + stage.cost;
                obj.objective = obj.objective + stage.objective;
                prev_fe = stage.stage(end);
            end
        end

        function createComplementarityConstraints(obj)
            import casadi.*
            model = obj.model;
            settings = obj.settings;
            dims = obj.dims;

            sigma_p = obj.sigma_p;
            s_elastic = obj.s_elastic;

            psi_fun = settings.psi_fun;

            if settings.elasticity_mode == ElasticityMode.NONE
                sigma = sigma_p;
            else
                sigma = s_elastic;
            end

            g_cross_comp = SX([]);
            % TODO Implement other modes
            if ~settings.use_fesd || settings.cross_comp_mode < 11
                % Do nothing, handled at the FE or stage level
                return
            elseif settings.cross_comp_mode == 11
                for r=1:dims.n_sys
                    g_r = 0;
                    nz_r = [];
                    for stage=obj.stages
                        for fe=stage.stage
                            pairs = fe.cross_comp_pairs(:, :, r);
                            expr_cell = cellfun(@(pair) apply_psi(pair, psi_fun, sigma), pairs, 'uni', false);
                            nonzeros = cellfun(@(x) vector_is_zero(x), expr_cell, 'uni', 0);
                            if size(vertcat(expr_cell{:}), 1) == 0
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
                            idx = exprs.sparsity().find();
                            if numel(idx) == 0
                                idx = [];
                            end
                            g_r = g_r + exprs(idx);
                            nz_r = nz_r + nonzeros;
                        end
                    end
                    g_r = scale_sigma(g_r, sigma, nz_r);
                    g_cross_comp = vertcat(g_cross_comp, g_r);
                end
            elseif settings.cross_comp_mode == 12
                for r=1:dims.n_sys
                    g_r = 0;
                    nz_r = [];
                    for stage=obj.stages
                        for fe=stage.stage
                            pairs = fe.cross_comp_pairs(:, :, r);
                            expr_cell = cellfun(@(pair) apply_psi(pair, psi_fun, sigma), pairs, 'uni', false);
                            nonzeros = cellfun(@(x) vector_is_zero(x), expr_cell, 'uni', 0);
                            if size(vertcat(expr_cell{:}), 1) == 0
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
                            if isempty(nz_r)
                                nz_r = zeros(size(nonzeros));
                            end
                            idx = exprs.sparsity().find();
                            if numel(idx) == 0
                                idx = [];
                            end
                            g_r = g_r + exprs(idx);
                            nz_r = nz_r + nonzeros;
                        end
                    end
                    g_r = scale_sigma(g_r, sigma, nz_r);
                    g_cross_comp = vertcat(g_cross_comp, g_r);
                end
            end

            % If We need to add a cost from the reformulation do that as needed;
            if settings.mpcc_mode == MpccMode.ell_1_penalty % this implies bilinear
                cost = 0;
                for r=1:dims.n_sys
                    for stage=obj.stages
                        for fe=stage.stage
                            pairs = fe.cross_comp_pairs(:, :, r);
                            expr_cell = cellfun(@(pair) apply_psi(pair, @(a,b,t) a.*b, 0), pairs, 'uni', false);
                            expr = sum1(sum2([expr_cell{:}]));
                            cost = cost + expr;
                        end
                    end
                end
                if settings.objective_scaling_direct
                    obj.cost = obj.cost + (1/sigma_p)*cost;
                else
                    obj.cost = sigma_p*obj.cost + cost;
                end
            else
                g_comp = g_cross_comp;
                n_comp = length(g_cross_comp);

                [g_comp_lb, g_comp_ub, g_comp] = generate_mpcc_relaxation_bounds(g_comp, settings);

                % Add reformulated constraints
                obj.addConstraint(g_comp, g_comp_lb, g_comp_ub);

            end
        end

        % TODO this should be private
        function addControlStage(obj, stage)
            w_len = length(obj.w);

            obj.addPrimalVector(stage.w, stage.lbw, stage.ubw, stage.w0);

            obj.ind_h = [obj.ind_h, increment_indices(stage.ind_h,w_len)];
            obj.ind_u = [obj.ind_u, {stage.ind_u+w_len}];
            obj.ind_sot = [obj.ind_sot, stage.ind_sot+w_len];
            obj.ind_x(stage.ctrl_idx, :, :) = increment_indices(stage.ind_x, w_len);
            obj.ind_v(stage.ctrl_idx, :, :) = increment_indices(stage.ind_v, w_len);
            obj.ind_theta(stage.ctrl_idx, :, :) = increment_indices(stage.ind_theta, w_len);
            obj.ind_lam(stage.ctrl_idx, :, :) = increment_indices(stage.ind_lam, w_len);
            obj.ind_mu(stage.ctrl_idx, :, :) = increment_indices(stage.ind_mu, w_len);
            obj.ind_alpha(stage.ctrl_idx, :, :) = increment_indices(stage.ind_alpha, w_len);
            obj.ind_lambda_n(stage.ctrl_idx, :, :) = increment_indices(stage.ind_lambda_n, w_len);
            obj.ind_lambda_p(stage.ctrl_idx, :, :) = increment_indices(stage.ind_lambda_p, w_len);
            obj.ind_beta(stage.ctrl_idx, :, :) = increment_indices(stage.ind_beta, w_len);
            obj.ind_theta_step(stage.ctrl_idx, :, :) = increment_indices(stage.ind_theta_step, w_len);
            obj.ind_z(stage.ctrl_idx, :, :) = increment_indices(stage.ind_z, w_len);
            obj.ind_nu_lift = [obj.ind_nu_lift, increment_indices(stage.ind_nu_lift, w_len)];

            obj.ind_lambda_normal(stage.ctrl_idx, :, :) = increment_indices(stage.ind_lambda_normal, w_len);
            obj.ind_lambda_tangent(stage.ctrl_idx, :, :) = increment_indices(stage.ind_lambda_tangent, w_len);
            obj.ind_y_gap(stage.ctrl_idx, :, :) = increment_indices(stage.ind_y_gap, w_len);
            obj.ind_gamma(stage.ctrl_idx, :, :) = increment_indices(stage.ind_gamma, w_len);
            obj.ind_beta_conic(stage.ctrl_idx, :, :) = increment_indices(stage.ind_beta_conic, w_len);
            obj.ind_gamma_d(stage.ctrl_idx, :, :) = increment_indices(stage.ind_gamma_d, w_len);
            obj.ind_beta_d(stage.ctrl_idx, :, :) = increment_indices(stage.ind_beta_d, w_len);
            obj.ind_delta_d(stage.ctrl_idx, :, :) = increment_indices(stage.ind_delta_d, w_len);
            obj.ind_p_vt(stage.ctrl_idx, :, :) = increment_indices(stage.ind_p_vt, w_len);
            obj.ind_n_vt(stage.ctrl_idx, :, :) = increment_indices(stage.ind_n_vt, w_len);
            obj.ind_alpha_vt(stage.ctrl_idx, :, :) = increment_indices(stage.ind_alpha_vt, w_len);

            obj.ind_x_left_bp(stage.ctrl_idx, :) = increment_indices(stage.ind_x_left_bp, w_len);
            obj.ind_Y_gap(stage.ctrl_idx, :) = increment_indices(stage.ind_Y_gap, w_len);
            obj.ind_Lambda_normal(stage.ctrl_idx, :) = increment_indices(stage.ind_Lambda_normal, w_len);
            obj.ind_Lambda_tangent(stage.ctrl_idx, :) = increment_indices(stage.ind_Lambda_tangent, w_len);
            obj.ind_Gamma(stage.ctrl_idx, :) = increment_indices(stage.ind_Gamma, w_len);
            obj.ind_Beta_conic(stage.ctrl_idx, :) = increment_indices(stage.ind_Beta_conic, w_len);
            obj.ind_Gamma_d(stage.ctrl_idx, :) = increment_indices(stage.ind_Gamma_d, w_len);
            obj.ind_Beta_d(stage.ctrl_idx, :) = increment_indices(stage.ind_Beta_d, w_len);
            obj.ind_Delta_d(stage.ctrl_idx, :) = increment_indices(stage.ind_Delta_d, w_len);
            obj.ind_L_vn(stage.ctrl_idx, :) = increment_indices(stage.ind_L_vn, w_len);
            obj.ind_P_vt(stage.ctrl_idx, :) = increment_indices(stage.ind_P_vt, w_len);
            obj.ind_N_vt(stage.ctrl_idx, :) = increment_indices(stage.ind_N_vt, w_len);
            obj.ind_Alpha_vt(stage.ctrl_idx, :) = increment_indices(stage.ind_Alpha_vt, w_len);

            obj.addConstraint(stage.g, stage.lbg, stage.ubg);
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

        function cc_vector = get.cc_vector(obj)
            cc_vector = [];
            for stage=obj.stages
                for fe=stage.stage
                    cc_vector = vertcat(cc_vector, fe.all_comp_pairs(:, 1) .* fe.all_comp_pairs(:, 2));
                end
            end
        end

        function u = get.u(obj)
            u = cellfun(@(u) obj.w(u), obj.ind_u, 'UniformOutput', false);
        end

        function sot = get.sot(obj)
            sot = cellfun(@(sot) obj.w(sot), obj.ind_sot, 'UniformOutput', false);
        end

        function nu_vector = get.nu_vector(obj)
            nu_vector = [];
            for k=obj.settings.N_stages
                stage = obj.stages(k);
                for fe=stage.stage
                    nu_vector = vertcat(nu_vector,fe.nu_vector);
                end
            end
        end

        function ind_z_all = get.ind_z_all(obj)
            ind_z_all = [flatten_ind(obj.ind_theta(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_lam(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_mu(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_alpha(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_lambda_n(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_lambda_p(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_beta(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_theta_step(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_z(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_lambda_normal(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_lambda_tangent(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_y_gap(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_gamma(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_beta_conic(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_gamma_d(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_beta_d(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_delta_d(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_p_vt(:,:,1:obj.dims.n_s))
                         flatten_ind(obj.ind_p_vt(:,:,1:obj.dims.n_s))
                        ];
            ind_z_all = sort(ind_z_all);
        end

        function ind_x_all = get.ind_x_all(obj)
            ind_x_all = [obj.ind_x0.'; flatten_ind(obj.ind_x)];
        end

        function print(obj,filename)
            if exist('filename')
                delete(filename);
                fileID = fopen(filename, 'w');
            else
                fileID = 1;
            end
            fprintf(fileID, "i\tlbg\t\t ubg\t\t g_expr\n");
            for i = 1:length(obj.lbg)
                expr_str = formattedDisplayText(obj.g(i));
                fprintf(fileID, "%d\t%.2e\t%.2e\t%s\n", i, obj.lbg(i), obj.ubg(i), expr_str);
            end

            fprintf(fileID, "\nw\t\t\tw0\t\tlbw\t\tubw\n");
            for i = 1:length(obj.lbw)
                % keyboard
                expr_str = pad(formattedDisplayText(obj.w(i)), 20);
                lb_str = pad(sprintf('%.2e', obj.lbw(i)), 10);
                fprintf(fileID, "%s\t%.2e\t%s\t%.2e\t\n", expr_str, obj.w0(i), lb_str, obj.ubw(i));
            end

            fprintf(fileID, '\nobjective\n');
            fprintf(fileID, strcat(formattedDisplayText(obj.cost), '\n'));
        end
    end
end
