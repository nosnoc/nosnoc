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
        ind_z
        % Stewart
        ind_theta
        ind_lam
        ind_mu
        % Step
        ind_alpha
        ind_lambda_n
        ind_lambda_p
        ind_beta
        ind_theta_step
        % CLS
        ind_x_left_bp
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
        ind_Y_gap
        ind_Lambda_normal
        ind_Lambda_tangent
        %
        ind_Gamma
        ind_Beta_conic
        %
        ind_Gamma_d
        ind_Beta_d
        ind_Delta_d
        %
        ind_L_vn
        ind_P_vt
        ind_N_vt
        ind_Alpha_vt
        % misc
        ind_h
        ind_nu_lift
        ind_elastic
        ind_comp_lift

        n_cont
        n_discont
        n_indep

        ind_eq
        ind_ineq
        ind_comp

        cross_comp_pairs
        all_comp_pairs
        n_comp_components
        ind_std_comp      % NOTE: EXPERIMENTAL, doesn't work for CLS, or c_n != 1 integration schemes

        ctrl_idx
        fe_idx

        model
        problem_options
        dims

        u

        T_final

        prev_fe
    end

    properties(Dependent, SetAccess=private, Hidden)
        x
        v
        % Stewart
        lam
        mu
        % Step
        alpha
        lambda_n
        lambda_p
        % CLS
        lambda_normal
        lambda_tangent
        y_gap
        % conic
        gamma
        beta_conic
        % poly
        gamma_d
        beta_d
        delta_d
        % variables related to conic
        p_vt
        n_vt
        alpha_vt
        % variables only at element boundary
        Y_Gap
        Lambda_normal
        Lambda_tangent
        Gamma
        Gamma_d
        Beta_conic
        Delta_d
        P_vn
        N_vn
        P_vt
        N_vt
        Alpha_vt
        % misc
        nu_lift
        h

        elastic

        nu_vector
    end

    properties(SetAccess=private)
        cross_comp_discont_0 % cross comp variables that are discont (e.g. theta in stewart)
    end

    methods
        function obj = FiniteElement(prev_fe, problem_options, model, dims, ctrl_idx, fe_idx, varargin)
            import casadi.*
            obj@NosnocFormulationObject();

            p = inputParser();
            p.FunctionName = 'FiniteElement';

            % TODO: add checks.
            addRequired(p, 'prev_fe');
            addRequired(p, 'problem_options');
            addRequired(p, 'model');
            addRequired(p, 'dims');
            addRequired(p, 'ctrl_idx');
            addRequired(p, 'fe_idx');
            addOptional(p, 'T_final',[]);
            parse(p, prev_fe, problem_options, model, dims, ctrl_idx, fe_idx, varargin{:});

            % store inputs
            obj.ctrl_idx = ctrl_idx;
            obj.fe_idx = fe_idx;

            obj.problem_options = problem_options;
            obj.model = model;
            obj.dims = dims;

            obj.prev_fe = prev_fe;

            % TODO Combine this with the below if
            if problem_options.right_boundary_point_explicit || problem_options.dcs_mode == 'CLS'
                rbp_allowance = 0;
            else
                rbp_allowance = 1;
            end

            if problem_options.dcs_mode == "CLS" && ~problem_options.right_boundary_point_explicit
                rbp_x_only = 1;
                obj.ind_x = cell(dims.n_s+1, 1);
            else
                rbp_x_only = 0;
                obj.ind_x = cell(dims.n_s+rbp_allowance, 1);
            end

            right_ygap = 0;
            if problem_options.dcs_mode == "CLS" && ~problem_options.right_boundary_point_explicit
                right_ygap = 1;
            end

            obj.ind_v = cell(dims.n_s, 1);
            obj.ind_z = cell(dims.n_s+rbp_allowance, 1);
            % Stewart
            obj.ind_theta = cell(dims.n_s+rbp_allowance, dims.n_sys);
            obj.ind_lam = cell(dims.n_s+rbp_allowance, dims.n_sys);
            obj.ind_mu = cell(dims.n_s+rbp_allowance, dims.n_sys);
            % Step
            obj.ind_alpha = cell(dims.n_s+rbp_allowance, dims.n_sys);
            obj.ind_lambda_n = cell(dims.n_s+rbp_allowance, dims.n_sys);
            obj.ind_lambda_p = cell(dims.n_s+rbp_allowance, dims.n_sys);
            obj.ind_theta_step = cell(dims.n_s+rbp_allowance, 1);
            obj.ind_beta = cell(dims.n_s+rbp_allowance, 1);
            % CLS
            obj.ind_lambda_normal = cell(dims.n_s+rbp_allowance,dims.n_sys);
            obj.ind_lambda_tangent = cell(dims.n_s+rbp_allowance,dims.n_sys);

            obj.ind_y_gap = cell(dims.n_s+right_ygap+rbp_allowance,dims.n_sys);

            % friction multipliers and lifting
            % conic
            obj.ind_gamma = cell(dims.n_s+rbp_allowance,dims.n_sys);
            obj.ind_beta_conic = cell(dims.n_s+rbp_allowance,dims.n_sys);
            % poly
            obj.ind_gamma_d = cell(dims.n_s+rbp_allowance,dims.n_sys);
            obj.ind_beta_d = cell(dims.n_s+rbp_allowance,dims.n_sys);
            obj.ind_delta_d = cell(dims.n_s+rbp_allowance,dims.n_sys);
            % variables related to conic
            obj.ind_p_vt = cell(dims.n_s+rbp_allowance,dims.n_sys);
            obj.ind_n_vt = cell(dims.n_s+rbp_allowance,dims.n_sys);
            obj.ind_alpha_vt = cell(dims.n_s+rbp_allowance,dims.n_sys);
            % variables only at element boundary
            obj.ind_Lambda_normal = cell(1,1);
            obj.ind_Y_gap = cell(1,1);
            obj.ind_Lambda_tangent = cell(1,1);
            obj.ind_Gamma = cell(1,1);
            obj.ind_Beta_d = cell(1,1);
            obj.ind_Gamma_d = cell(1,1);
            obj.ind_Beta_conic = cell(1,1);
            obj.ind_Delta_d = cell(1,1);
            obj.ind_L_vn = cell(1,1);
            obj.ind_P_vt = cell(1,1);
            obj.ind_N_vt = cell(1,1);
            obj.ind_Alpha_vt = cell(1,1);
            obj.ind_x_left_bp = cell(1,1);

            % misc
            obj.ind_h = [];
            obj.ind_elastic = [];
            
            if problem_options.use_fesd
                h = define_casadi_symbolic(problem_options.casadi_symbolic_mode, ['h_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1)]);
                h_ctrl_stage = problem_options.T/problem_options.N_stages;
                h0 = h_ctrl_stage / problem_options.N_finite_elements(ctrl_idx);
                ubh = (1 + problem_options.gamma_h) * h0;
                lbh = (1 - problem_options.gamma_h) * h0;
                if problem_options.time_rescaling && ~problem_options.use_speed_of_time_variables
                    % if only time_rescaling is true, speed of time and step size all lumped together, e.g., \hat{h}_{k,i} = s_n * h_{k,i}, hence the bounds need to be extended.
                    ubh = (1+problem_options.gamma_h)*h0*problem_options.s_sot_max;
                    lbh = (1-problem_options.gamma_h)*h0/problem_options.s_sot_min;
                end
                obj.addVariable(h, 'h', lbh, ubh, h0);
            end
            obj.T_final = p.Results.T_final;

            % Left boundary point needed for dcs_mode = cls (corresponding to t^+ (post impacts))
            if problem_options.dcs_mode == "CLS"
                x = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                    ['X_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_0'],...
                    dims.n_x);
                obj.addVariable(x,...
                    'x_left_bp',...
                    model.lbx,...
                    model.ubx,...
                    model.x0,...
                    1);
            end
            % RK stage stuff
            for ii = 1:dims.n_s
                % state / state derivative variables
                % TODO @Anton better v initialization.
                if (problem_options.rk_representation == RKRepresentation.differential ||...
                        problem_options.rk_representation == RKRepresentation.differential_lift_x)
                    v = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                        ['V_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)'],...
                        dims.n_x);
                    obj.addVariable(v,...
                        'v',...
                        -inf * ones(dims.n_x,1),...
                        inf * ones(dims.n_x,1),...
                        zeros(dims.n_x,1),...
                        ii);
                end
                if (problem_options.rk_representation == RKRepresentation.integral ||...
                        problem_options.rk_representation == RKRepresentation.differential_lift_x)
                    if problem_options.x_box_at_stg
                        lbx = model.lbx;
                        ubx = model.ubx;
                    elseif problem_options.x_box_at_fe && ii == dims.n_s && problem_options.right_boundary_point_explicit
                        lbx = model.lbx;
                        ubx = model.ubx;
                    elseif fe_idx == problem_options.N_finite_elements(ctrl_idx) && ii == dims.n_s && problem_options.right_boundary_point_explicit
                        lbx = model.lbx;
                        ubx = model.ubx;
                    else
                        lbx = -inf * ones(dims.n_x, 1);
                        ubx = inf * ones(dims.n_x, 1);
                    end
                    x = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
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
                % TODO @Anton revive lp_based_guess here
                if problem_options.dcs_mode == DcsMode.Stewart
                    % add thetas
                    for ij = 1:dims.n_sys
                        initial_theta = (1/dims.n_f_sys(ij)) * ones(dims.n_f_sys(ij),1);
                        theta = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                            ['theta_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            dims.n_f_sys(ij));
                        obj.addVariable(theta,...
                            'theta',...
                            zeros(dims.n_f_sys(ij), 1),...
                            inf * ones(dims.n_f_sys(ij), 1),...
                            initial_theta,...
                            ii,...
                            ij);
                    end
                    % add lambdas
                    for ij = 1:dims.n_sys
                        initial_lambda = ones(dims.n_f_sys(ij), 1);
                        lam = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                            ['lambda_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            dims.n_f_sys(ij));
                        obj.addVariable(lam,...
                            'lam',...
                            zeros(dims.n_f_sys(ij), 1),...
                            inf * ones(dims.n_f_sys(ij), 1),...
                            initial_lambda,...
                            ii,...
                            ij);
                    end
                    % add mu
                    for ij = 1:dims.n_sys
                        mu = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                            ['mu_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            1);
                        obj.addVariable(mu,...
                            'mu',...
                            -inf,...
                            inf,...
                            0.5,...
                            ii,...
                            ij);
                    end
                elseif problem_options.dcs_mode == DcsMode.Heaviside
                    % add alpha
                    for ij = 1:dims.n_sys
                        alpha = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                            ['alpha_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            dims.n_c_sys(ij));
                        obj.addVariable(alpha,...
                            'alpha',...
                            zeros(dims.n_c_sys(ij), 1),...
                            ones(dims.n_c_sys(ij), 1),...
                            problem_options.initial_alpha * ones(dims.n_c_sys(ij), 1),...
                            ii,...
                            ij);
                    end
                    % add lambda_n and lambda_p
                    for ij = 1:dims.n_sys
                        lambda_n = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                            ['lambda_n_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            dims.n_c_sys(ij));
                        lambda_p = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                            ['lambda_p_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii) '_' num2str(ij)],...
                            dims.n_c_sys(ij));
                        obj.addVariable(lambda_n,...
                            'lambda_n',...
                            zeros(dims.n_c_sys(ij), 1),...
                            inf * ones(dims.n_c_sys(ij), 1),...
                            0.5 * ones(dims.n_c_sys(ij), 1),...
                            ii,...
                            ij);
                        obj.addVariable(lambda_p,...
                            'lambda_p',...
                            zeros(dims.n_c_sys(ij), 1),...
                            inf * ones(dims.n_c_sys(ij), 1),...
                            0.5 * ones(dims.n_c_sys(ij), 1),...
                            ii,...
                            ij);
                    end
                    % TODO: Clean this up (maybe as a function to reduce the indent level.
                    if problem_options.time_freezing_inelastic
                        %                         if ~problem_options.pss_lift_step_functions
                        theta_step = define_casadi_symbolic(problem_options.casadi_symbolic_mode, ['theta_step_'  num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],dims.n_theta_step);
                        beta = define_casadi_symbolic(problem_options.casadi_symbolic_mode, ['beta_'  num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],dims.n_beta);
                        beta_guess = full(model.g_lift_beta_fun(problem_options.initial_alpha*ones(dims.n_alpha,1)));
                        theta_step_guess = full(model.g_lift_theta_step_fun(problem_options.initial_alpha*ones(dims.n_alpha,1),beta_guess));
                        obj.addVariable(theta_step,...
                            'theta_step',...
                            -inf*ones(dims.n_theta_step,1),...
                            inf*ones(dims.n_theta_step,1),...
                            theta_step_guess,...
                            ii);
                        if problem_options.pss_lift_step_functions
                            obj.addVariable(beta,...
                                'beta',...
                                -inf*ones(dims.n_beta,1),...
                                inf*ones(dims.n_beta,1),...
                                beta_guess,...
                                ii);
                        end
                    end
                elseif problem_options.dcs_mode == DcsMode.CLS
                    lambda_normal = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                        ['lambda_normal_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                        dims.n_contacts);
                    obj.addVariable(lambda_normal,...
                        'lambda_normal',...
                        zeros(dims.n_contacts,1),...
                        inf * ones(dims.n_contacts, 1),...
                        0*ones(dims.n_contacts, 1),...
                        ii,1);
                    y_gap = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                        ['y_gap_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                        dims.n_contacts);
                    obj.addVariable(y_gap,...
                        'y_gap',...
                        zeros(dims.n_contacts,1),...
                        inf * ones(dims.n_contacts, 1),...
                        0*ones(dims.n_contacts, 1),...
                        ii,1);
                    if model.friction_exists
                        lambda_tangent = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                            ['lambda_tangent_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                            dims.n_tangents);

                        if isequal(problem_options.friction_model,'Polyhedral')
                            obj.addVariable(lambda_tangent,...
                                'lambda_tangent',...
                                zeros(dims.n_tangents,1),...
                                inf * ones(dims.n_tangents, 1),...
                                ones(dims.n_tangents, 1),...
                                ii,1);

                            gamma_d = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                                ['gamma_d_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                dims.n_contacts);
                            obj.addVariable(gamma_d,...
                                'gamma_d',...
                                zeros(dims.n_contacts,1),...
                                inf * ones(dims.n_contacts, 1),...
                                ones(dims.n_contacts, 1),...
                                ii,1);

                            beta_d = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                                ['beta_d_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                dims.n_contacts);
                            obj.addVariable(beta_d,...
                                'beta_d',...
                                zeros(dims.n_contacts,1),...
                                inf * ones(dims.n_contacts, 1),...
                                ones(dims.n_contacts, 1),...
                                ii,1);

                            delta_d = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                                ['delta_d_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                dims.n_tangents);

                            obj.addVariable(delta_d,...
                                'delta_d',...
                                zeros(dims.n_tangents,1),...
                                inf * ones(dims.n_tangents, 1),...
                                ones(dims.n_tangents, 1),...
                                ii,1);
                        end
                        if isequal(problem_options.friction_model,'Conic')
                            obj.addVariable(lambda_tangent,...
                                'lambda_tangent',...
                                -inf*ones(dims.n_tangents,1),...
                                inf * ones(dims.n_tangents, 1),...
                                ones(dims.n_tangents, 1),...
                                ii,1);

                            gamma = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                                ['gamma_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                dims.n_contacts);

                            obj.addVariable(gamma,...
                                'gamma',...
                                zeros(dims.n_contacts,1),...
                                inf * ones(dims.n_contacts, 1),...
                                ones(dims.n_contacts, 1),...
                                ii,1);

                            beta = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                                ['beta_conic_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                dims.n_contacts);

                            obj.addVariable(beta,...
                                'beta_conic',...
                                zeros(dims.n_contacts,1),...
                                inf * ones(dims.n_contacts, 1),...
                                ones(dims.n_contacts, 1),...
                                ii,1);

                            switch problem_options.conic_model_switch_handling
                                case 'Plain'
                                    % no extra vars
                                case 'Abs'
                                    p_vt = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                                        ['p_vt_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                        dims.n_tangents);
                                    obj.addVariable(p_vt,...
                                        'p_vt',...
                                        zeros(dims.n_tangents,1),...
                                        inf * ones(dims.n_tangents, 1),...
                                        ones(dims.n_tangents, 1),...
                                        ii,1);
                                    n_vt = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                                        ['n_vt_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                        dims.n_tangents);
                                    obj.addVariable(n_vt,...
                                        'n_vt',...
                                        zeros(dims.n_tangents,1),...
                                        inf * ones(dims.n_tangents, 1),...
                                        ones(dims.n_tangents, 1),...
                                        ii,1);
                                case 'Lp'
                                    p_vt = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                                        ['p_vt_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                        dims.n_tangents);
                                    obj.addVariable(p_vt,...
                                        'p_vt',...
                                        zeros(dims.n_tangents,1),...
                                        inf * ones(dims.n_tangents, 1),...
                                        ones(dims.n_tangents, 1),...
                                        ii,1);
                                    n_vt = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                                        ['n_vt_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                        dims.n_tangents);
                                    obj.addVariable(n_vt,...
                                        'n_vt',...
                                        zeros(dims.n_tangents,1),...
                                        inf * ones(dims.n_tangents, 1),...
                                        ones(dims.n_tangents, 1),...
                                        ii,1);
                                    alpha_vt = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                                        ['alpha_vt_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                                        dims.n_tangents);
                                    obj.addVariable(alpha_vt,...
                                        'alpha_vt',...
                                        zeros(dims.n_tangents,1),...
                                        inf * ones(dims.n_tangents, 1),...
                                        ones(dims.n_tangents, 1),...
                                        ii,1);
                            end
                        end
                    end
                end
                % add user variables
                z = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                    ['z_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                    dims.n_z);
                obj.addVariable(z,...
                    'z',...
                    model.lbz,...
                    model.ubz,...
                    model.z0,...
                    ii);
            end

            if problem_options.dcs_mode == DcsMode.CLS && (obj.fe_idx ~= 1 || ~problem_options.no_initial_impacts)
                %  IMPULSE VARIABLES
                Lambda_normal = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                    ['Lambda_normal_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1)],...
                    dims.n_contacts);
                obj.addVariable(Lambda_normal,...
                    'Lambda_normal',...
                    zeros(dims.n_contacts,1),...
                    inf * ones(dims.n_contacts, 1),...
                    0*ones(dims.n_contacts, 1),1);

                L_vn = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                    ['L_vn' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                    dims.n_contacts);
                obj.addVariable(L_vn,...
                    'L_vn',...
                    -inf*ones(dims.n_contacts,1),...
                    inf * ones(dims.n_contacts, 1),...
                    ones(dims.n_contacts, 1),1);

                Y_gap = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                    ['Y_gap_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                    dims.n_contacts);
                obj.addVariable(Y_gap,...
                    'Y_gap',...
                    zeros(dims.n_contacts,1),...
                    inf * ones(dims.n_contacts, 1),...
                    0*ones(dims.n_contacts, 1),1);
                if model.friction_exists
                    Lambda_tangent = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                        ['Lambda_tangent' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                        dims.n_tangents);

                    if problem_options.friction_model == FrictionModel.Polyhedral
                        obj.addVariable(Lambda_tangent,...
                            'Lambda_tangent',...
                            zeros(dims.n_tangents,1),...
                            inf * ones(dims.n_tangents, 1),...
                            ones(dims.n_tangents, 1),1);

                        Gamma_d = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                            ['Gamma_d' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                            dims.n_contacts);
                        obj.addVariable(Gamma_d,...
                            'Gamma_d',...
                            zeros(dims.n_contacts,1),...
                            inf * ones(dims.n_contacts, 1),...
                            ones(dims.n_contacts, 1),1);

                        Beta_d = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                            ['Beta_d' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                            dims.n_contacts);
                        obj.addVariable(Beta_d,...
                            'Beta_d',...
                            zeros(dims.n_contacts,1),...
                            inf * ones(dims.n_contacts, 1),...
                            ones(dims.n_contacts, 1),1);

                        Delta_d = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                            ['Delta_d' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                            dims.n_tangents);

                        obj.addVariable(Delta_d,...
                            'Delta_d',...
                            zeros(dims.n_tangents,1),...
                            inf * ones(dims.n_tangents, 1),...
                            ones(dims.n_tangents, 1),1);
                    end
                    if problem_options.friction_model == FrictionModel.Conic
                        obj.addVariable(Lambda_tangent,...
                            'Lambda_tangent',...
                            -inf*ones(dims.n_tangents,1),...
                            inf * ones(dims.n_tangents, 1),...
                            ones(dims.n_tangents, 1),1);

                        Gamma = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                            ['Gamma' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                            dims.n_contacts);

                        obj.addVariable(Gamma,...
                            'Gamma',...
                            zeros(dims.n_contacts,1),...
                            inf * ones(dims.n_contacts, 1),...
                            ones(dims.n_contacts, 1),1);

                        Beta_conic = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                            ['Beta_conic' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                            dims.n_contacts);

                        obj.addVariable(Beta_conic,...
                            'Beta_conic',...
                            zeros(dims.n_contacts,1),...
                            inf * ones(dims.n_contacts, 1),...
                            ones(dims.n_contacts, 1),1);

                        switch problem_options.conic_model_switch_handling
                            case 'Plain'
                                % no extra vars
                            case 'Abs'
                                P_vt = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                                    ['P_vt' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                                    dims.n_tangents);
                                obj.addVariable(P_vt,...
                                    'P_vt',...
                                    zeros(dims.n_tangents,1),...
                                    inf * ones(dims.n_tangents, 1),...
                                    ones(dims.n_tangents, 1),1);
                                N_vt = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                                    ['N_vt' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                                    dims.n_tangents);
                                obj.addVariable(N_vt,...
                                    'N_vt',...
                                    zeros(dims.n_tangents,1),...
                                    inf * ones(dims.n_tangents, 1),...
                                    ones(dims.n_tangents, 1),1);
                            case 'Lp'
                                P_vt = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                                    ['P_vt' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                                    dims.n_tangents);
                                obj.addVariable(P_vt,...
                                    'P_vt',...
                                    zeros(dims.n_tangents,1),...
                                    inf * ones(dims.n_tangents, 1),...
                                    ones(dims.n_tangents, 1),1);
                                N_vt = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                                    ['N_vt' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                                    dims.n_tangents);
                                obj.addVariable(N_vt,...
                                    'N_vt',...
                                    zeros(dims.n_tangents,1),...
                                    inf * ones(dims.n_tangents, 1),...
                                    ones(dims.n_tangents, 1),1);
                                Alpha_vt = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                                    ['Alpha_vt' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) ],...
                                    dims.n_tangents);
                                obj.addVariable(Alpha_vt,...
                                    'Alpha_vt',...
                                    zeros(dims.n_tangents,1),...
                                    inf * ones(dims.n_tangents, 1),...
                                    ones(dims.n_tangents, 1),1);
                        end
                    end
                end
            end


            % Add right boundary points if needed
            if ~problem_options.right_boundary_point_explicit
                if problem_options.dcs_mode == DcsMode.Stewart
                    % add lambdas
                    for ij = 1:dims.n_sys
                        lam = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                            ['lambda_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_end_' num2str(ij)],...
                            dims.n_f_sys(ij));
                        obj.addVariable(lam,...
                            'lam',...
                            zeros(dims.n_f_sys(ij), 1),...
                            inf * ones(dims.n_f_sys(ij), 1),...
                            ones(dims.n_f_sys(ij), 1),...
                            dims.n_s+1,...
                            ij);
                    end

                    % add mu
                    for ij = 1:dims.n_sys
                        mu = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                            ['mu_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_end_' num2str(ij)],...
                            1);
                        obj.addVariable(mu,...
                            'mu',...
                            -inf * ones(1),...
                            inf * ones(1),...
                            0.5 * ones(1),...
                            dims.n_s+1,...
                            ij);
                    end

                elseif problem_options.dcs_mode == DcsMode.Heaviside
                    % add lambda_n
                    for ij = 1:dims.n_sys
                        lambda_n = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                            ['lambda_n_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_end_' num2str(ij)],...
                            dims.n_c_sys(ij));
                        obj.addVariable(lambda_n,...
                            'lambda_n',...
                            zeros(dims.n_c_sys(ij), 1),...
                            inf * ones(dims.n_c_sys(ij), 1),...
                            0.5 * ones(dims.n_c_sys(ij), 1),...
                            dims.n_s+1,...
                            ij);
                    end

                    % add lambda_p
                    for ij = 1:dims.n_sys
                        lambda_p = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                            ['lambda_p_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_end_' num2str(ij)],...
                            dims.n_c_sys(ij));
                        obj.addVariable(lambda_p,...
                            'lambda_p',...
                            zeros(dims.n_c_sys(ij), 1),...
                            inf * ones(dims.n_c_sys(ij), 1),...
                            0.5 * ones(dims.n_c_sys(ij), 1),...
                            dims.n_s+1,...
                            ij);
                    end
                elseif problem_options.dcs_mode == DcsMode.CLS && obj.fe_idx == problem_options.N_finite_elements(obj.ctrl_idx)
                    % only for last finite element (in ctrl stage).
                    y_gap_end = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                        ['y_gap_end_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1) '_' num2str(ii)],...
                        dims.n_contacts);
                    obj.addVariable(y_gap_end,...
                        'y_gap',...
                        zeros(dims.n_contacts,1),...
                        inf * ones(dims.n_contacts, 1),...
                        0*ones(dims.n_contacts, 1),...
                        dims.n_s+1, 1);
                end
            end

            if (~problem_options.right_boundary_point_explicit ||...
                    problem_options.rk_representation == RKRepresentation.differential)
                if problem_options.x_box_at_stg || problem_options.x_box_at_fe || fe_idx == problem_options.N_finite_elements(ctrl_idx)
                    lbx = model.lbx;
                    ubx = model.ubx;
                else
                    lbx = -inf * ones(dims.n_x);
                    ubx = inf * ones(dims.n_x);
                end

                % add final X variables
                x = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                    ['X_end_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1)],...
                    dims.n_x);
                obj.addVariable(x,...
                    'x',...
                    lbx,...
                    ubx,...
                    model.x0,...
                    dims.n_s+rbp_allowance+rbp_x_only);
            end
            if strcmpi(problem_options.step_equilibration, 'direct_homotopy_lift')
                nu_lift = define_casadi_symbolic(problem_options.casadi_symbolic_mode,...
                    ['nu_lift_' num2str(ctrl_idx-1) '_' num2str(fe_idx-1)],...
                    1);
                obj.addVariable(nu_lift,...
                                'nu_lift',...
                                -inf,...
                                inf,...
                                1);
            end
        end
        
        function cross_comp_pairs = getCrossCompPairs(obj)
            import casadi.*
            model = obj.model;
            problem_options = obj.problem_options;
            dims = obj.dims;
            prev_fe = obj.prev_fe;
            
            switch problem_options.dcs_mode
                case DcsMode.Stewart
                    rbp_allowance = ~problem_options.right_boundary_point_explicit;
                    n_cont = dims.n_s+1+rbp_allowance;
                    n_discont = dims.n_s;
                    n_indep = dims.n_sys;
                    cross_comp_pairs = cell(n_cont, n_discont, n_indep);
                    if problem_options.use_fesd
                        % Generate cross complementarity pairs from within the finite element.
                        for j=2:n_cont
                            for jj=1:n_discont
                                for r=1:n_indep
                                    cross_comp_pairs{j,jj,r} = [obj.w(obj.ind_lam{j-1,r}), obj.w(obj.ind_theta{jj,r})];
                                end
                            end
                        end

                        % Generate cross complementarity pairs with the end of the last finite element.
                        for jj=1:n_discont
                            for r=1:n_indep
                                cross_comp_pairs{1,jj,r} = [prev_fe.w(prev_fe.ind_lam{end,r}), obj.w(obj.ind_theta{jj,r})];
                            end
                        end
                    else
                        for j=1:n_discont
                            for r=1:n_indep
                                cross_comp_pairs{j,j,r} = [obj.w(obj.ind_lam{j,r}), obj.w(obj.ind_theta{j,r})];
                            end
                        end
                    end
                case DcsMode.Heaviside
                    rbp_allowance = ~problem_options.right_boundary_point_explicit;
                    n_cont = dims.n_s+1+rbp_allowance;
                    n_discont = dims.n_s;
                    n_indep = dims.n_sys;
                    cross_comp_pairs = cell(n_cont, n_discont, n_indep);
                    if problem_options.use_fesd
                        % Generate cross complementarity pairs from within the finite element.
                        for j=2:n_cont
                            for jj=1:n_discont
                                for r=1:n_indep
                                    cross_comp_pairs{j,jj,r} = [vertcat(obj.w(obj.ind_lambda_n{j-1,r}),obj.w(obj.ind_lambda_p{j-1,r})),...
                                        vertcat(obj.w(obj.ind_alpha{jj,r}), ones(dims.n_c_sys(r),1)-obj.w(obj.ind_alpha{jj,r}))];
                                end
                            end
                        end

                        % Generate cross complementarity pairs with the end of the last finite element.
                        for jj=1:n_discont
                            for r=1:n_indep
                                cross_comp_pairs{1,jj,r} = [vertcat(prev_fe.w(prev_fe.ind_lambda_n{end,r}),prev_fe.w(prev_fe.ind_lambda_p{end,r})),...
                                    vertcat(obj.w(obj.ind_alpha{jj,r}), ones(dims.n_c_sys(r),1)-obj.w(obj.ind_alpha{jj,r}))];
                            end
                        end
                    else 
                        for j=1:n_discont
                            for r=1:n_indep
                                cross_comp_pairs{j,j,r} = [vertcat(obj.w(obj.ind_lambda_n{j,r}),obj.w(obj.ind_lambda_p{j,r})),...
                                    vertcat(obj.w(obj.ind_alpha{j,r}), ones(dims.n_c_sys(r),1)-obj.w(obj.ind_alpha{j,r}))];
                            end
                        end
                    end
                case DcsMode.CLS
                    n_cont = dims.n_s+2;
                    n_discont = dims.n_s+2;
                    n_indep = 1;
                    cross_comp_pairs = cell(n_cont, n_discont, n_indep);
                    % Generate cross complementarity pairs with the end of the last finite element.
                    % NOTE: we start from the 3rd column here because we include the complementarity pairs with the impulse
                    %       variables in the 2nd column and the complementarity pairs with the previous finite element in
                    %       the 1st column. This is slightly different than in Step or Stewart where we only have the previous
                    %       finite element to worry about. We iterate over the values for n_cont and n_discont therefore we
                    %       must subtract the buffer added for these when we index into the index sets.
                    if problem_options.use_fesd
                        for j=3:n_cont
                            for jj=3:n_discont
                                pairs = [];

                                pairs = vertcat(pairs, [obj.w(obj.ind_y_gap{j-2,1}), obj.w(obj.ind_lambda_normal{jj-2,1})]);
                                if model.friction_exists
                                    if problem_options.friction_model == FrictionModel.Conic
                                        pairs = vertcat(pairs, [obj.w(obj.ind_gamma{j-2,1}), obj.w(obj.ind_beta_conic{jj-2,1})]);
                                        switch problem_options.conic_model_switch_handling
                                            case 'Plain'
                                                % no extra expr
                                            case 'Abs'
                                                pairs = vertcat(pairs, [obj.w(obj.ind_p_vt{j-2,1}), obj.w(obj.ind_n_vt{jj-2,1})]);
                                            case 'Lp'
                                                pairs = vertcat(pairs, [obj.w(obj.ind_p_vt{j-2,1}), obj.w(obj.ind_alpha_vt{jj-2,1})]);
                                                pairs = vertcat(pairs, [obj.w(obj.ind_n_vt{j-2,1}), ones(dims.n_tangents,1)-obj.w(obj.ind_alpha_vt{jj-2,1})]);
                                        end
                                    elseif problem_options.friction_model == FrictionModel.Polyhedral
                                        pairs = vertcat(pairs, [obj.w(obj.ind_delta_d{j-2,1}), obj.w(obj.ind_lambda_tangent{jj-2,1})]);
                                        pairs = vertcat(pairs, [obj.w(obj.ind_gamma_d{j-2,1}), obj.w(obj.ind_beta_d{jj-2,1})]);
                                    end
                                end
                                cross_comp_pairs{j,jj,1} = pairs;
                            end
                        end

                        % Here we generate the cross complemetarities with the previous FE and the impulse variables.
                        for jj=1:n_discont-2
                            % Ignore the previous fe lambda normal if it is not the end point.
                            if problem_options.right_boundary_point_explicit || (obj.fe_idx == 1 && obj.ctrl_idx == 1)
                                cross_comp_pairs{1,jj+2} = [prev_fe.w(prev_fe.ind_y_gap{end,1}), obj.w(obj.ind_lambda_normal{jj,1})];
                            end
                            if model.friction_exists
                                if problem_options.friction_model == FrictionModel.Conic
                                    if problem_options.right_boundary_point_explicit || (obj.fe_idx == 1 && obj.ctrl_idx == 1)
                                        cross_comp_pairs{1,jj+2} = vertcat(cross_comp_pairs{1,jj+2}, [prev_fe.w(prev_fe.ind_gamma{end,1}), obj.w(obj.ind_beta_conic{jj,1})]);
                                    end
                                    switch problem_options.conic_model_switch_handling
                                        case 'Plain'
                                            % no extra expr
                                        case 'Abs'
                                            cross_comp_pairs{1,jj+2} = vertcat(cross_comp_pairs{1,jj+2}, [prev_fe.w(prev_fe.ind_p_vt{end,1}), obj.w(obj.ind_n_vt{jj,1})]);
                                        case 'Lp'
                                            cross_comp_pairs{1,jj+2} = vertcat(cross_comp_pairs{1,jj+2}, [prev_fe.w(prev_fe.ind_p_vt{end,1}), obj.w(obj.ind_alpha_vt{jj,1})]);
                                            cross_comp_pairs{1,jj+2} = vertcat(cross_comp_pairs{1,jj+2}, [prev_fe.w(prev_fe.ind_n_vt{end,1}), ones(dims.n_tangents,1)-obj.w(obj.ind_alpha_vt{jj,1})]);
                                    end
                                elseif problem_options.friction_model == FrictionModel.Polyhedral
                                    cross_comp_pairs{1,jj+2} = vertcat(cross_comp_pairs{1,jj+2}, [prev_fe.w(prev_fe.ind_delta_d{end,1}), obj.w(obj.ind_lambda_tangent{jj,1})]);
                                    cross_comp_pairs{1,jj+2} = vertcat(cross_comp_pairs{1,jj+2}, [prev_fe.w(prev_fe.ind_gamma_d{end,1}), obj.w(obj.ind_beta_d{jj,1})]);
                                end
                            end
                            if (obj.fe_idx ~= 1 || ~problem_options.no_initial_impacts)
                                comp = [obj.w(obj.ind_Y_gap{1}), obj.w(obj.ind_lambda_normal{jj,1})]; % TODO maybe we need the other cross comps?
                                size_comp = size(comp, 1);
                                size_full_comp = size(cross_comp_pairs{end,jj+2});
                                sparsity = Sparsity.rowcol([0:(size_comp-1)],[0,1],size_full_comp(1),size_full_comp(2));
                                cross_comp_pairs{2,jj+2} = define_casadi_symbolic(problem_options.casadi_symbolic_mode, 'dummy', size_full_comp, sparsity);
                                cross_comp_pairs{2,jj+2}(1:size_comp, :) = comp;
                            end
                        end
                    else % no fesd
                        for jj=3:n_discont
                            pairs = [];

                            pairs = vertcat(pairs, [obj.w(obj.ind_y_gap{jj-2,1}), obj.w(obj.ind_lambda_normal{jj-2,1})]);
                            if model.friction_exists
                                if problem_options.friction_model == FrictionModel.Conic
                                    pairs = vertcat(pairs, [obj.w(obj.ind_gamma{jj-2,1}), obj.w(obj.ind_beta_conic{jj-2,1})]);
                                    switch problem_options.conic_model_switch_handling
                                        case 'Plain'
                                            % no extra expr
                                        case 'Abs'
                                            pairs = vertcat(pairs, [obj.w(obj.ind_p_vt{jj-2,1}), obj.w(obj.ind_n_vt{jj-2,1})]);
                                        case 'Lp'
                                            pairs = vertcat(pairs, [obj.w(obj.ind_p_vt{jj-2,1}), obj.w(obj.ind_alpha_vt{jj-2,1})]);
                                            pairs = vertcat(pairs, [obj.w(obj.ind_n_vt{jj-2,1}), ones(dims.n_tangents,1)-obj.w(obj.ind_alpha_vt{jj-2,1})]);
                                    end
                                elseif problem_options.friction_model == FrictionModel.Polyhedral
                                    pairs = vertcat(pairs, [obj.w(obj.ind_delta_d{jj-2,1}), obj.w(obj.ind_lambda_tangent{jj-2,1})]);
                                    pairs = vertcat(pairs, [obj.w(obj.ind_gamma_d{jj-2,1}), obj.w(obj.ind_beta_d{jj-2,1})]);
                                end
                            end
                            cross_comp_pairs{jj,jj,1} = pairs;
                        end
                    end
            end
            obj.n_cont = n_cont;
            obj.n_discont = n_discont;
            obj.n_indep = n_indep;
            obj.cross_comp_pairs = cross_comp_pairs;
        end

        function h = get.h(obj)
            if obj.problem_options.use_fesd
                h = obj.w(obj.ind_h);
            elseif obj.problem_options.time_optimal_problem && ~obj.problem_options.use_speed_of_time_variables
                h = obj.T_final/(obj.problem_options.N_stages*obj.problem_options.N_finite_elements(obj.ctrl_idx));
            else
                h = obj.problem_options.T/(obj.problem_options.N_stages*obj.problem_options.N_finite_elements(obj.ctrl_idx));
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

            % TODO handle independent subsystems?
            if obj.problem_options.use_fesd && obj.fe_idx > 1
                switch obj.problem_options.dcs_mode
                    case DcsMode.Stewart
                        lam_F = cellfun(@(x) obj.w(x), obj.ind_lam, 'uni', false);
                        lam_B = cellfun(@(x) obj.prev_fe.w(x), obj.prev_fe.ind_lam, 'uni', false);
                        theta_F = cellfun(@(x) obj.w(x), obj.ind_theta, 'uni', false);
                        theta_B = cellfun(@(x) obj.prev_fe.w(x), obj.prev_fe.ind_theta, 'uni', false);
                        sigma_cont_F = sum2(horzcat(lam_F{:}));
                        sigma_cont_B = sum2(horzcat(lam_B{:}));
                        sigma_discont_F = sum2(horzcat(theta_F{:}));
                        sigma_discont_B = sum2(horzcat(theta_B{:}));

                        pi_cont = sigma_cont_B .* sigma_cont_F;
                        pi_discont = sigma_discont_B .* sigma_discont_F;
                        nu = pi_cont + pi_discont;
                    case DcsMode.Heaviside
                        lam_p_F = cellfun(@(x) obj.w(x), obj.ind_lambda_p, 'uni', false);
                        lam_p_B = cellfun(@(x) obj.prev_fe.w(x), obj.prev_fe.ind_lambda_p, 'uni', false);
                        lam_n_F = cellfun(@(x) obj.w(x), obj.ind_lambda_n, 'uni', false);
                        lam_n_B = cellfun(@(x) obj.prev_fe.w(x), obj.prev_fe.ind_lambda_n, 'uni', false);
                        alpha_F = cellfun(@(x) obj.w(x), obj.ind_alpha, 'uni', false);
                        alpha_B = cellfun(@(x) obj.prev_fe.w(x), obj.prev_fe.ind_alpha, 'uni', false);

                        sigma_cont_F = vertcat(sum2(horzcat(lam_n_F{:})), sum2(horzcat(lam_p_F{:})));
                        sigma_cont_B = vertcat(sum2(horzcat(lam_n_B{:})), sum2(horzcat(lam_p_B{:})));
                        sigma_discont_F = vertcat(sum2(horzcat(alpha_F{:})), sum2(1-horzcat(alpha_F{:})));
                        sigma_discont_B = vertcat(sum2(horzcat(alpha_B{:})), sum2(1-horzcat(alpha_B{:})));

                        pi_cont = sigma_cont_B .* sigma_cont_F;
                        pi_discont = sigma_discont_B .* sigma_discont_F;
                        nu = pi_cont + pi_discont;
                    case DcsMode.CLS
                        f_c_F = cellfun(@(x) obj.w(x), obj.ind_y_gap, 'uni', false);
                        f_c_B = cellfun(@(x) obj.prev_fe.w(x), obj.prev_fe.ind_y_gap, 'uni', false);
                        l_n_F = cellfun(@(x) obj.w(x), obj.ind_lambda_normal, 'uni', false);
                        l_n_B = cellfun(@(x) obj.prev_fe.w(x), obj.prev_fe.ind_lambda_normal, 'uni', false);
                        
                        sigma_f_c_F = obj.w(obj.ind_Y_gap{1})+sum2(horzcat(f_c_F{:}));
                        % If we don't have an old Y gap due to no initial impacts ignore it in backward sum
                        if isempty(obj.prev_fe.ind_Y_gap{1})
                            sigma_f_c_B = sum2(horzcat(f_c_B{:}));
                        else
                            sigma_f_c_B = obj.prev_fe.w(obj.prev_fe.ind_Y_gap{1})+sum2(horzcat(f_c_B{:}));
                        end
                        sigma_l_n_F = sum2(horzcat(l_n_F{:}));
                        sigma_l_n_B = sum2(horzcat(l_n_B{:}));

                        pi_f_c = obj.w(obj.ind_Y_gap{1}).*sigma_f_c_F.*sigma_f_c_B;
                        pi_l_n = sigma_l_n_F .* sigma_l_n_B;
                        kappa = pi_f_c + pi_l_n;
                        if obj.model.friction_exists
                            % Go one by one for each contact (TODO vectorize)
                            zeta = SX.ones(obj.model.dims.n_contacts, 1);
                            for ii=1:obj.model.dims.n_contacts
                                if obj.model.mu_f(ii)
                                    if strcmp(obj.problem_options.friction_model, 'Conic')
                                        xi_p_F = cellfun(@(x) obj.w(x), obj.ind_p_vt, 'uni', false);
                                        xi_p_B = cellfun(@(x) obj.prev_fe.w(x), obj.prev_fe.ind_p_vt, 'uni', false);
                                        xi_n_F = cellfun(@(x) obj.w(x), obj.ind_n_vt, 'uni', false);
                                        xi_n_B = cellfun(@(x) obj.prev_fe.w(x), obj.prev_fe.ind_n_vt, 'uni', false);
                                        beta_F = cellfun(@(x) obj.w(x), obj.ind_beta_conic, 'uni', false);
                                        beta_B = cellfun(@(x) obj.prev_fe.w(x), obj.prev_fe.ind_beta_conic, 'uni', false);

                                        sigma_xi_p_F = obj.w(obj.ind_P_vt{1})+sum2(horzcat(xi_p_F{:}));
                                        sigma_xi_n_F = obj.w(obj.ind_N_vt{1})+sum2(horzcat(xi_n_F{:}));
                                        if isempty(obj.prev_fe.ind_P_vt{1})
                                            sigma_xi_p_B = sum2(horzcat(xi_p_B{:}));
                                            sigma_xi_n_B = sum2(horzcat(xi_n_B{:}));
                                        else
                                            sigma_xi_p_B = obj.prev_fe.w(obj.prev_fe.ind_P_vt{1})+sum2(horzcat(xi_p_B{:}));
                                            sigma_xi_n_B = obj.prev_fe.w(obj.prev_fe.ind_N_vt{1})+sum2(horzcat(xi_n_B{:}));
                                        end
                                        sigma_beta_F = sum2(horzcat(beta_F{:}));
                                        sigma_beta_B = sum2(horzcat(beta_B{:}));

                                        pi_xi_p = sigma_xi_p_F .* sigma_xi_p_B;
                                        pi_xi_n = sigma_xi_n_F .* sigma_xi_n_B;
                                        pi_beta = sigma_beta_F .* sigma_beta_B;

                                        zeta(ii) = (sigma_f_c_B + sigma_f_c_F) + (sum(pi_xi_p + pi_xi_n) + pi_beta);
                                    else
                                        % TODO: Verify
                                        gamma_F = cellfun(@(x) obj.w(x), obj.ind_gamma_d, 'uni', false);
                                        gamma_B = cellfun(@(x) obj.prev_fe.w(x), obj.prev_fe.ind_gamma_d, 'uni', false);
                                        beta_F = cellfun(@(x) obj.w(x), obj.ind_beta_d, 'uni', false);
                                        beta_B = cellfun(@(x) obj.prev_fe.w(x), obj.prev_fe.ind_beta_d, 'uni', false);

                                        sigma_gamma_F = sum2(horzcat(gamma_F{:}));
                                        sigma_gamma_B = sum2(horzcat(gamma_B{:}));
                                        sigma_beta_F = sum2(horzcat(beta_F{:}));
                                        sigma_beta_B = sum2(horzcat(beta_B{:}));
                                        
                                        pi_gamma = sigma_gamma_F.*sigma_gamma_B;
                                        pi_beta = sigma_beta_F.*sigma_beta_B;

                                        zeta(ii) = 1;%(sigma_f_c_B + sigma_f_c_F) + (sum(pi_gamma) + sum(pi_beta));
                                    end
                                end
                            end
                            nu = kappa.*zeta;
                        else
                            nu = kappa;
                        end
                end
                
                nu_vector = 1;
                for jjj=1:length(nu)
                    nu_vector = nu_vector * nu(jjj);
                end
            else
                nu_vector = [];
            end
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
                   [obj.ind_z{stage}],...
                   [obj.ind_lambda_normal{stage}],...
                   [obj.ind_y_gap{stage}],...
                   [obj.ind_lambda_tangent{stage}],...
                   [obj.ind_gamma{stage}],...
                   [obj.ind_beta_conic{stage}],...
                   [obj.ind_gamma_d{stage}],...
                   [obj.ind_beta_d{stage}],...
                   [obj.ind_delta_d{stage}],...
                   [obj.ind_p_vt{stage}],...
                   [obj.ind_n_vt{stage}],...
                   [obj.ind_alpha_vt{stage}],...
                  ];

            z = obj.w(idx);
        end

        function forwardSimulation(obj, ocp, Uk, s_sot, p_stage)
            model = obj.model;
            problem_options = obj.problem_options;
            dims = obj.dims;

            obj.u = Uk;
            % left bondary point
            if problem_options.dcs_mode == "CLS"
                % do continuity on x and impulse equtions
                X_k0 = obj.w(obj.ind_x_left_bp{1}); % corresponds to post impact t^+
                X_k = obj.prev_fe.x{end}; % corresponds to pre impact t^-
                Q_k0  = X_k0(1:dims.n_q);
                V_k0  = X_k0(dims.n_q+1:end);
                Q_k  = X_k(1:dims.n_q);
                V_k  = X_k(dims.n_q+1:end);
                % junction equations
                obj.addConstraint(Q_k0-Q_k);
                %  Z_impulse_k = [obj.w(obj.ind_Lambda_normal{1}); obj.w(obj.ind_Y_gap{1});obj.w(obj.ind_P_vn{1});obj.w(obj.ind_N_vn{1});...
                %  obj.w(obj.ind_Lambda_tangent{1}); obj.w(obj.ind_Gamma_d{1}); obj.w(obj.ind_Beta_d{1}); obj.w(obj.ind_Delta_d{1}); ...
                %  obj.w(obj.ind_Gamma{1}); obj.w(obj.ind_Beta_conic{1}); obj.w(obj.ind_P_vt{1}); obj.w(obj.ind_N_vt{1}); obj.w(obj.ind_Alpha_vt{1})];
                Z_impulse_k = [obj.w(obj.ind_Lambda_normal{1}); obj.w(obj.ind_Y_gap{1});obj.w(obj.ind_L_vn{1});...
                               obj.w(obj.ind_Lambda_tangent{1}); obj.w(obj.ind_Gamma_d{1}); obj.w(obj.ind_Beta_d{1}); obj.w(obj.ind_Delta_d{1}); ...
                    obj.w(obj.ind_Gamma{1}); obj.w(obj.ind_Beta_conic{1}); obj.w(obj.ind_P_vt{1}); obj.w(obj.ind_N_vt{1}); obj.w(obj.ind_Alpha_vt{1})];
                if (obj.fe_idx ~= 1 || ~problem_options.no_initial_impacts)
                    obj.addConstraint(model.g_impulse_fun(Q_k0,V_k0,V_k,Z_impulse_k));
                else
                    obj.addConstraint(V_k0-V_k);
                end
                % additional y_eps constraint
                if problem_options.eps_cls > 0
                    x_eps = vertcat(Q_k0 + obj.h * problem_options.eps_cls * V_k0, V_k0);
                    obj.addConstraint( model.f_c_fun(x_eps), 0, inf);
                end
            else
                X_k0 = obj.prev_fe.x{end};
            end

            if problem_options.rk_representation == RKRepresentation.integral
                X_ki = obj.x;
                Xk_end = problem_options.D_irk(1) * X_k0;
            elseif problem_options.rk_representation == RKRepresentation.differential
                X_ki = {};
                for j = 1:dims.n_s
                    x_temp = X_k0;
                    for r = 1:dims.n_s
                        x_temp = x_temp + obj.h*problem_options.A_irk(j,r)*obj.v{r};
                    end
                    X_ki = [X_ki {x_temp}];
                end
                X_ki = [X_ki, {obj.x{end}}];
                Xk_end = X_k0;
            elseif problem_options.rk_representation == RKRepresentation.differential_lift_x
                X_ki = obj.x;
                Xk_end = X_k0;
                for j = 1:dims.n_s
                    x_temp = X_k0;
                    for r = 1:dims.n_s
                        x_temp = x_temp + obj.h*problem_options.A_irk(j,r)*obj.v{r};
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
                if problem_options.rk_representation == RKRepresentation.integral
                    xj = problem_options.C_irk(1, j+1) * X_k0;
                    for r = 1:dims.n_s
                        xj = xj + problem_options.C_irk(r+1, j+1) * X_ki{r};
                    end
                    Xk_end = Xk_end + problem_options.D_irk(j+1) * X_ki{j};
                    obj.addConstraint(obj.h * fj - xj);
                    if problem_options.cost_integration
                        obj.augmented_objective = obj.augmented_objective + problem_options.B_irk(j+1) * obj.h * qj;
                        obj.objective = obj.objective + problem_options.B_irk(j+1) * obj.h * qj;
                    end
                else
                    Xk_end = Xk_end + obj.h * problem_options.b_irk(j) * obj.v{j};
                    obj.addConstraint(fj - obj.v{j});
                    if problem_options.cost_integration
                         obj.augmented_objective = obj.augmented_objective + problem_options.b_irk(j) * obj.h * qj;
                         obj.objective = obj.objective + problem_options.b_irk(j) * obj.h * qj;
                    end
                end
            end
            if ~problem_options.cost_integration && problem_options.N_finite_elements(obj.ctrl_idx) == obj.fe_idx
                if problem_options.right_boundary_point_explicit
                    [~, q] = model.f_x_fun(X_ki{dims.n_s}, obj.rkStageZ(dims.n_s), Uk, p_stage, model.v_global);
                else
                    [~, q] = model.f_x_fun(X_ki{dims.n_s+1}, obj.rkStageZ(dims.n_s+1), Uk, p_stage, model.v_global);
                end
                obj.augmented_objective = obj.augmented_objective + q;
                obj.objective = obj.objective + q;
            end

            % nonlinear inequality.
            % TODO: do this cleaner
            if (~isempty(model.g_path) &&...
                (obj.fe_idx == problem_options.N_finite_elements(obj.ctrl_idx) || problem_options.g_path_at_fe))
                endidx = dims.n_s + ~problem_options.right_boundary_point_explicit;
                obj.addConstraint(model.g_path_fun(obj.x{end},obj.rkStageZ(endidx),Uk,p_stage,model.v_global), model.g_path_lb, model.g_path_ub);
            end
            for j=1:dims.n_s-problem_options.right_boundary_point_explicit
                % TODO: there has to be a better way to do this.
                if ~isempty(model.g_path) && problem_options.g_path_at_stg
                    obj.addConstraint(model.g_path_fun(X_ki{j},obj.rkStageZ(j),Uk,p_stage,model.v_global), model.g_path_lb, model.g_path_ub);
                end
            end

            % end constraints
            if (~problem_options.right_boundary_point_explicit ||...
                    problem_options.rk_representation == RKRepresentation.differential)
                obj.addConstraint(Xk_end - obj.x{end});
            end
            if (~problem_options.right_boundary_point_explicit &&...
                problem_options.use_fesd &&...
                problem_options.dcs_mode ~= DcsMode.CLS)

                % TODO verify this.
                obj.addConstraint(model.g_switching_fun(obj.x{end}, obj.rkStageZ(dims.n_s+1), Uk, p_stage));
            end
            % y_gap_end
            if (~problem_options.right_boundary_point_explicit &&...
                problem_options.N_finite_elements(obj.ctrl_idx) == obj.fe_idx &&...
                problem_options.dcs_mode == DcsMode.CLS)
                obj.addConstraint(model.f_c_fun(obj.x{end}) - obj.w(obj.ind_y_gap{end}))
            end
        end

        function createComplementarityConstraints(obj, p_stage)
            import casadi.*
            model = obj.model;
            problem_options = obj.problem_options;
            dims = obj.dims;

            % TODO: Lift all pairs that are not scalar.
            % TODO: how to determine pairs that are not scalar?
            g_path_comp_pairs = [];

            % path complementarities
            if (~isempty(model.g_comp_path) &&...
                (obj.fe_idx == problem_options.N_finite_elements(obj.ctrl_idx) || problem_options.g_path_at_fe))
                endidx = dims.n_s + ~problem_options.right_boundary_point_explicit;
                pairs = model.g_comp_path_fun(obj.prev_fe.x{end}, obj.rkStageZ(endidx), obj.u, p_stage, model.v_global);
                g_path_comp_pairs = vertcat(g_path_comp_pairs, pairs);
            end
            for j=1:dims.n_s-problem_options.right_boundary_point_explicit
                % TODO: there has to be a better way to do this.
                if ~isempty(model.g_comp_path) && problem_options.g_path_at_stg
                    % FIXME this actually doesn't work for pure differential mode.
                    pairs = model.g_comp_path_fun(obj.x{j}, obj.rkStageZ(j), obj.u, p_stage, model.v_global);
                    g_path_comp_pairs = vertcat(g_path_comp_pairs, pairs);
                end
            end
            
            
            impulse_pairs = [];
            if problem_options.dcs_mode == DcsMode.CLS && (obj.fe_idx ~= 1 || ~problem_options.no_initial_impacts)
                % comp condts
                % Y_gap comp. to Lambda_Normal+P_vn+N_vn;..
                % P_vn comp.to N_vn;
                %  if conic:
                    % Gamma comp. to Beta
                    % if abs: P_vt comp.to N_vt;
                    % if lp:  P_vt comp.to e - Alpha_vt, N_vt comp.to Alpha_vt
                %  if polyhedral
                   % Delta comp. to Lambda_tangent
                   % Gamma comp. to Beta;
                
                Y_gap = obj.w(obj.ind_Y_gap{1});
                Lambda_normal = obj.w(obj.ind_Lambda_normal{1});
                % P_vn = obj.w(obj.ind_P_vn{1});
                % N_vn = obj.w(obj.ind_N_vn{1});
                L_vn = obj.w(obj.ind_L_vn{1});
                P_vt = obj.w(obj.ind_P_vt{1});
                N_vt = obj.w(obj.ind_N_vt{1});
                Gamma = obj.w(obj.ind_Gamma{1});
                Beta_conic = obj.w(obj.ind_Beta_conic{1});
                Lambda_tangent = obj.w(obj.ind_Lambda_tangent{1});
                Beta_d = obj.w(obj.ind_Beta_d{1});
                Alpha_vt = obj.w(obj.ind_Alpha_vt{1});
                Delta_d = obj.w(obj.ind_Delta_d{1});
                Gamma_d = obj.w(obj.ind_Gamma_d{1});
                
                % impulse_pairs = vertcat(impulse_pairs, [Lambda_normal, (Y_gap+P_vn+N_vn)]);
                % impulse_pairs = vertcat(impulse_pairs, [P_vn, N_vn]);
                % constraint on position and velocity level separated
                impulse_pairs = vertcat(impulse_pairs, [Lambda_normal, (Y_gap)]);
                impulse_pairs = vertcat(impulse_pairs, [Lambda_normal, L_vn]);
                impulse_pairs = vertcat(impulse_pairs, [Lambda_normal, -L_vn]);
                if model.friction_exists
                    if problem_options.friction_model == FrictionModel.Conic
                        impulse_pairs = vertcat(impulse_pairs, [Gamma,Beta_conic]);
                        switch problem_options.conic_model_switch_handling
                          case ConicModelSwitchHandling.Plain
                            % no extra comps
                          case ConicModelSwitchHandling.Abs
                            impulse_pairs = vertcat(impulse_pairs, [P_vt,N_vt]);
                          case ConicModelSwitchHandling.Lp
                            impulse_pairs = vertcat(impulse_pairs, [P_vt,1-Alpha_vt]);
                            impulse_pairs = vertcat(impulse_pairs, [Alpha_vt,N_vt]);
                        end
                    elseif problem_options.friction_model == FrictionModel.Polyhedral
                        impulse_pairs = vertcat(impulse_pairs, [Delta_d,Lambda_tangent]);
                        impulse_pairs = vertcat(impulse_pairs, [Gamma_d,Beta_d]);
                    end
                end
                % if we are not in the first element we need to also cross comp Ygap  with lambdas of prev fe
                % TODO this should be done cleaner but for now this is fine.
                if (obj.fe_idx ~= 1) && ~problem_options.right_boundary_point_explicit
                    for ii=1:obj.prev_fe.n_discont-2
                        impulse_pairs = vertcat(impulse_pairs, [Y_gap, obj.prev_fe.w(obj.prev_fe.ind_lambda_normal{ii, 1})]);
                    end
                end

                if (obj.fe_idx == problem_options.N_finite_elements(obj.ctrl_idx)) && ~problem_options.right_boundary_point_explicit
                    for ii=1:problem_options.n_s
                        impulse_pairs = vertcat(impulse_pairs, [obj.w(obj.ind_y_gap{end}), obj.w(obj.ind_lambda_normal{ii, 1})]);
                    end
                end
            end
            
            cross_comp_pairs = obj.getCrossCompPairs();

            cross_comp_aggregated = [];
            ind_std_comp = [];

            sigma_scale = 1; % TODO scale properly
                             % apply psi
            g_cross_comp = [];
            lbg_cross_comp = [];
            ubg_cross_comp = [];
            if problem_options.cross_comp_mode == CrossCompMode.STAGE_STAGE || ~problem_options.use_fesd
                % NOTE: this is a hack.
                bool_cells = cellfun(@(x) false(size(x,1),1), cross_comp_pairs, 'uni', false);
                if problem_options.dcs_mode == "CLS"
                    for ii=3:obj.n_discont
                        bool_cells{ii, ii} = ~bool_cells{ii, ii};
                    end
                else
                    for ii=1:obj.n_discont
                        for r=1:obj.n_indep
                            bool_cells{1+ii, ii,r} = ~bool_cells{1+ii, ii,r};
                        end

                    end
                end
                obj.ind_std_comp = vertcat(bool_cells{:});
                cross_comp_aggregated = vertcat(cross_comp_pairs{:});
            elseif problem_options.cross_comp_mode == CrossCompMode.FE_STAGE
                a = [];
                b = [];
                for j=1:obj.n_cont
                    for r=1:obj.n_indep
                        pairs = cross_comp_pairs(j, :, r);
                        n_pair = size(pairs{end}, 1);
                        if all(cellfun(@isempty, pairs))
                            continue
                        end
                        idx = find(~cellfun(@isempty,pairs),1);
                        cont = pairs{idx}(:,1);
                        discont = cellfun(@(pair) vertcat(pair, zeros(n_pair - length(pair),2)),pairs, 'uni', false);
                        discont = cellfun(@(pair) pair(:,2), discont, 'uni', false);
                        b = [b;sum2([discont{:}])];
                        a = [a;cont];
                    end
                end
                obj.ind_std_comp = true(size(a)); % TODO This may be wrong
                cross_comp_aggregated = [a,b];
            elseif problem_options.cross_comp_mode == CrossCompMode.STAGE_FE
                a = [];
                b = [];
                for jj=1:obj.n_discont
                    for r=1:obj.n_indep
                        pairs = cross_comp_pairs(:, jj, r);
                        n_pair = size(pairs{end}, 1);
                        if all(cellfun(@isempty, pairs))
                            continue
                        end
                        idx = find(~cellfun(@isempty,pairs),1);
                        discont = pairs{idx}(:,2);
                        cont = cellfun(@(pair) vertcat(pair, zeros(n_pair - length(pair),2)),pairs, 'uni', false);
                        cont = cellfun(@(pair) pair(:,1), cont, 'uni', false);
                        b = [b;discont];
                        a = [a;sum2([cont{:}])];
                    end
                end
                obj.ind_std_comp = true(size(a)); % TODO This may be wrong
                cross_comp_aggregated = [a,b];
            elseif problem_options.cross_comp_mode == CrossCompMode.FE_FE
                a = [];
                b = [];
                for r=1:obj.n_indep
                    pairs = cross_comp_pairs(:, :, r);
                    if all(cellfun(@isempty, pairs))
                        continue
                    end
                    n_pair = size(pairs{end,end}, 1);
                    cont = cellfun(@(pair) vertcat(pair, zeros(n_pair - length(pair),2)),pairs(:,end), 'uni', false);
                    cont = cellfun(@(pair) pair(:,1), cont, 'uni', false);
                    discont = cellfun(@(pair) vertcat(pair, zeros(n_pair - length(pair),2)),pairs(end,:), 'uni', false);
                    discont = cellfun(@(pair) pair(:,2), discont, 'uni', false);
                    b = [b;sum2([discont{:}])];
                    a = [a;sum2([cont{:}])];
                end
                obj.ind_std_comp = true(size(a));
                cross_comp_aggregated = [a,b];
            end
            obj.ind_std_comp = vertcat(true(size(g_path_comp_pairs,1),1),true(size(impulse_pairs,1),1), obj.ind_std_comp);
            obj.all_comp_pairs = vertcat(g_path_comp_pairs, impulse_pairs, cross_comp_aggregated);

            if problem_options.lift_complementarities
                n_lift_comp = 0;
                ind_lift_comp = [];
                tmp = obj.all_comp_pairs(:);
                % figure out what we need to lift
                for ii=1:length(tmp)
                    if ~is_leaf(tmp(ii))
                        n_lift_comp = n_lift_comp+1;
                        ind_lift_comp = [ind_lift_comp, ii];
                    end
                end
                z_comp = define_casadi_symbolic(problem_options.casadi_symbolic_mode, ['z_comp_' num2str(obj.ctrl_idx-1) '_' num2str(obj.fe_idx-1)], n_lift_comp);
                if problem_options.lower_bound_comp_lift
                    lb = zeros(n_lift_comp,1);
                else
                    lb = -inf*ones(n_lift_comp,1);
                end
                obj.addVariable(z_comp,...
                    'comp_lift',...
                    lb,...
                    inf*ones(n_lift_comp, 1),...
                    ones(n_lift_comp, 1));
                if length(ind_lift_comp)
                    g_comp_lift = z_comp - tmp(ind_lift_comp);
                    obj.addConstraint(g_comp_lift);
                    obj.all_comp_pairs(ind_lift_comp) = z_comp;
                end
            end

            if problem_options.experimental_supervertical_form
                % lift G
                n_lift_comp = 0;
                ind_lift_comp = [];
                tmp = obj.all_comp_pairs(:,1);
                % figure out what we need to lift
                for ii=1:length(tmp)
                    n_lift_comp = n_lift_comp+1;
                    ind_lift_comp = [ind_lift_comp, ii];
                end
                z_comp = define_casadi_symbolic(problem_options.casadi_symbolic_mode, ['z_G_' num2str(obj.ctrl_idx-1) '_' num2str(obj.fe_idx-1)], n_lift_comp);
                if problem_options.lower_bound_comp_lift
                    lb = zeros(n_lift_comp,1);
                else
                    lb = -inf*ones(n_lift_comp,1);
                end
                obj.addVariable(z_comp,...
                    'comp_lift',...
                    lb,...
                    inf*ones(n_lift_comp, 1),...
                    ones(n_lift_comp, 1));
                if length(ind_lift_comp)
                    g_comp_lift = z_comp - tmp(ind_lift_comp);
                    obj.addConstraint(g_comp_lift);
                    obj.all_comp_pairs(ind_lift_comp,1) = z_comp;
                end

                % lift H
                n_lift_comp = 0;
                ind_lift_comp = [];
                tmp = obj.all_comp_pairs(:,2);
                % figure out what we need to lift
                for ii=1:length(tmp)
                    n_lift_comp = n_lift_comp+1;
                    ind_lift_comp = [ind_lift_comp, ii];
                end
                z_comp = define_casadi_symbolic(problem_options.casadi_symbolic_mode, ['z_H_' num2str(obj.ctrl_idx-1) '_' num2str(obj.fe_idx-1)], n_lift_comp);
                if problem_options.lower_bound_comp_lift
                    lb = zeros(n_lift_comp,1);
                else
                    lb = -inf*ones(n_lift_comp,1);
                end
                obj.addVariable(z_comp,...
                    'comp_lift',...
                    lb,...
                    inf*ones(n_lift_comp, 1),...
                    ones(n_lift_comp, 1));
                if length(ind_lift_comp)
                    g_comp_lift = z_comp - tmp(ind_lift_comp);
                    obj.addConstraint(g_comp_lift);
                    obj.all_comp_pairs(ind_lift_comp,2) = z_comp;
                end
            end
        end

        function stepEquilibration(obj, rho_h_p)
            import casadi.*
            model = obj.model;
            problem_options = obj.problem_options;
            dims = obj.dims;

            if ~problem_options.use_fesd
                return;
            end

            % only heuristic mean is done for first finite element
            if problem_options.step_equilibration == StepEquilibrationMode.heuristic_mean
                h_fe = problem_options.T / (sum(problem_options.N_finite_elements)); % TODO this may be a bad idea if using different N_fe. may want to issue warning in that case
                obj.augmented_objective = obj.augmented_objective + rho_h_p * (obj.h - h_fe).^2;
                return;
            elseif obj.fe_idx <= 1
                return;
            end

            % step equilibration modes that depend on previous FE.
            nu = obj.nu_vector;
            delta_h_ki = obj.h - obj.prev_fe.h;
            if problem_options.step_equilibration ==  StepEquilibrationMode.heuristic_diff
                obj.augmented_objective = obj.augmented_objective + rho_h_p * delta_h_ki.^2;
            elseif problem_options.step_equilibration == StepEquilibrationMode.l2_relaxed_scaled
                obj.augmented_objective = obj.augmented_objective + rho_h_p * tanh(nu/problem_options.step_equilibration_sigma) * delta_h_ki.^2;
            elseif problem_options.step_equilibration == StepEquilibrationMode.l2_relaxed
                obj.augmented_objective = obj.augmented_objective + rho_h_p * nu * delta_h_ki.^2;
            elseif problem_options.step_equilibration == StepEquilibrationMode.direct
                obj.addConstraint(nu*delta_h_ki, 0, 0);
                % TODO: how to do this if it is not part of the mpcc.
            elseif problem_options.step_equilibration == StepEquilibrationMode.direct_homotopy
                obj.addConstraint([nu*delta_h_ki-rho_h_p;-nu*delta_h_ki-rho_h_p],...
                    [-inf;-inf],...
                    [0;0]);
            elseif problem_options.step_equilibration == StepEquilibrationMode.direct_homotopy_lift
                obj.addConstraint([obj.nu_lift-nu;obj.nu_lift*delta_h_ki-rho_h_p;-obj.nu_lift*delta_h_ki-rho_h_p],...
                    [0;-inf;-inf],...
                    [0;0;0]);
            else
                error("Step equilibration mode not implemented");
            end
        end

        function json = jsonencode(obj, varargin)
            import casadi.*
            fe_struct = struct(obj);

            fe_struct = rmfield(fe_struct, 'prev_fe');
            fe_struct = rmfield(fe_struct, 'model');
            fe_struct = rmfield(fe_struct, 'dims');
            fe_struct = rmfield(fe_struct, 'problem_options');
            json = jsonencode(fe_struct);
        end
    end
end

