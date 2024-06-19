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
classdef FiniteElementZero < NosnocFormulationObject
    properties
        model
        settings
        % Index vectors
        ind_x
        % Stewart
        ind_lam
        % Step
        ind_lambda_n
        ind_lambda_p
        % CLS
        ind_y_gap
        % conic
        ind_gamma
        % poly
        ind_gamma_d
        ind_delta_d
        % variables related to conic
        ind_p_vt
        ind_n_vt
        % variables only at element boundary
        ind_Gamma
        ind_Gamma_d
        ind_Delta_d
%         ind_P_vn
%         ind_N_vn
        ind_L_vn

        % Parameter for initial conditions
        x0
    end

    properties(Dependent, SetAccess=private, Hidden)
        % Casadi symbolics for FE variables
        x
        cross_comp_cont_0
        cross_comp_cont_1
        cross_comp_cont_2
        %         lambda
        % CLS
        y_gap
        % conic
        gamma
        % poly
        gamma_d
        delta_d
        % variables related to conic
        p_vt
        n_vt
        %         variables only at element boundary
        %         TODO: for impulse at 0th
        Gamma
        Gamma_d
        Delta_d
        P_vn
        N_vn
    end

    methods
        function obj = FiniteElementZero(settings, dims, model)
            import casadi.*
            obj@NosnocFormulationObject();
            
            obj.model = model;
            obj.settings = settings;
            obj.ind_x = cell(1, 1);
            obj.ind_lam = cell(1,dims.n_sys);
            obj.ind_lambda_n = cell(1,dims.n_sys);
            obj.ind_lambda_p = cell(1,dims.n_sys);

            obj.ind_y_gap = cell(1,dims.n_sys);
            obj.ind_gamma = cell(1,dims.n_sys);
            obj.ind_gamma_d= cell(1,dims.n_sys);
            obj.ind_delta_d = cell(1,dims.n_sys);
            obj.ind_p_vt = cell(1,dims.n_sys);
            obj.ind_n_vt = cell(1,dims.n_sys);

            X0 = define_casadi_symbolic(settings.casadi_symbolic_mode, 'X0', dims.n_x); % variable
            obj.x0 = define_casadi_symbolic(settings.casadi_symbolic_mode, 'x0', dims.n_x); % Param

            if model.there_exist_free_x0
                for i=1:dims.n_x
                    if ~ismember(i, model.ind_free_x0)
                        obj.addConstraint(X0(i) - obj.x0(i),0,0);
                    end
                end
                x0_ub = inf*ones(dims.n_x,1);
                x0_lb = -inf*ones(dims.n_x,1);
                x0_ub(model.ind_free_x0) = ubx(model.ind_free_x0);
                x0_lb(model.ind_free_x0) = lbx(model.ind_free_x0);

                obj.addVariable(X0,...
                    'x',...
                    x0_lb,...
                    x0_ub,...
                    model.x0,...
                    1);
            else
                for i=1:dims.n_x
                    obj.addConstraint(X0(i) - obj.x0(i),0,0);
                end
                obj.addVariable(X0,...
                    'x',...
                    -inf*ones(dims.n_x, 1),...
                    inf*ones(dims.n_x, 1),...
                    model.x0,...
                    1);
            end

            % lambda00
            % TODO Use define_casadi_symbolic
            if settings.dcs_mode == DcsMode.Stewart
                for ij=1:dims.n_sys
                    lam0 = define_casadi_symbolic(settings.casadi_symbolic_mode, ['lambda00_' num2str(ij)], dims.n_f_sys(ij));
                    obj.addVariable(lam0,...
                        'lam',...
                        -inf * ones(dims.n_f_sys(ij),1),...
                        inf * ones(dims.n_f_sys(ij),1),...
                        ones(dims.n_f_sys(ij),1),...
                        1,...
                        ij);
                end
            elseif settings.dcs_mode == DcsMode.Heaviside
                for ij=1:dims.n_sys
                    lambda_n0 = define_casadi_symbolic(settings.casadi_symbolic_mode, ['lambda00_n_' num2str(ij)], dims.n_c_sys(ij));
                    lambda_p0 = define_casadi_symbolic(settings.casadi_symbolic_mode, ['lambda00_p_' num2str(ij)], dims.n_c_sys(ij));
                    obj.addVariable(lambda_n0,...
                        'lambda_n',...
                        zeros(dims.n_c_sys(ij),1),...
                        inf * ones(dims.n_c_sys(ij),1),...
                        0.5 * ones(dims.n_c_sys(ij),1),...
                        1,...
                        ij);
                    obj.addVariable(lambda_p0,...
                        'lambda_p',...
                        zeros(dims.n_c_sys(ij),1),...
                        inf * ones(dims.n_c_sys(ij),1),...
                        0.5 * ones(dims.n_c_sys(ij),1),...
                        1,...
                        ij);
                end

            elseif settings.dcs_mode == DcsMode.CLS
                % TODO: implement impulse equations for the zeroth element
                y_gap0 = define_casadi_symbolic(settings.casadi_symbolic_mode,'y_gap0', dims.n_contacts);
                obj.addVariable(y_gap0,...
                    'y_gap',...
                    zeros(dims.n_contacts,1),...
                    inf * ones(dims.n_contacts,1),...
                    ones(dims.n_contacts,1),...
                    1);
                if model.friction_exists
                    if isequal(settings.friction_model,'Polyhedral')

                        gamma_d0 = define_casadi_symbolic(settings.casadi_symbolic_mode,'gamma_d0', dims.n_contacts);
                        obj.addVariable(gamma_d0 ,...
                            'gamma_d',...
                            zeros(dims.n_contacts,1),...
                            inf * ones(dims.n_contacts,1),...
                            ones(dims.n_contacts,1),...
                            1);
                        delta_d0 = define_casadi_symbolic(settings.casadi_symbolic_mode,'delta_d0', dims.n_tangents);
                        obj.addVariable(delta_d0 ,...
                            'delta_d',...
                            zeros(dims.n_tangents,1),...
                            inf * ones(dims.n_tangents,1),...
                            ones(dims.n_tangents,1),...
                            1);
                        %
                    end
                    if isequal(settings.friction_model,'Conic')
                        gamma_0 = define_casadi_symbolic(settings.casadi_symbolic_mode,'gamma_0', dims.n_contacts);
                        obj.addVariable(gamma_0 ,...
                                        'gamma',...
                                        zeros(dims.n_contacts,1),...
                                        inf * ones(dims.n_contacts,1),...
                                        ones(dims.n_contacts,1),...
                                        1);
                        switch settings.conic_model_switch_handling
                          case 'Plain'
                            % no extra vars
                          case {'Abs','Lp'}
                            p_vt_0 = define_casadi_symbolic(settings.casadi_symbolic_mode,'p_vt_0', dims.n_tangents);
                            obj.addVariable(p_vt_0 ,...
                                            'p_vt',...
                                            zeros(dims.n_tangents,1),...
                                            inf * ones(dims.n_tangents,1),...
                                            ones(dims.n_tangents,1),...
                                            1);
                            n_vt_0 = define_casadi_symbolic(settings.casadi_symbolic_mode,'n_vt_0', dims.n_tangents);
                            obj.addVariable(n_vt_0 ,...
                                            'n_vt',...
                                            zeros(dims.n_tangents,1),...
                                            inf * ones(dims.n_tangents,1),...
                                            ones(dims.n_tangents,1),...
                                            1);
                        end
                    end
                end
            end

        end

        function cross_comp_cont_0 = get.cross_comp_cont_0(obj)
            import casadi.*
            %             grab = @(l, ln, lp) vertcat(obj.w(l), obj.w(ln), obj.w(lp));
            %             lambda = cellfun(grab, obj.ind_lam, obj.ind_lambda_n, obj.ind_lambda_p, 'UniformOutput', false);
            if obj.model.friction_exists && obj.settings.friction_model == 'Conic' && obj.settings.conic_model_switch_handling == 'Abs'
                grab = @(l, ln, lp ,yg, g, pvt, nvt, gd,dd) vertcat(obj.w(l), obj.w(ln), obj.w(lp), obj.w(yg), obj.w(g), obj.w(gd), obj.w(dd));
                cross_comp_cont_0 = cellfun(grab, obj.ind_lam, obj.ind_lambda_n, obj.ind_lambda_p,...
                    obj.ind_y_gap, obj.ind_gamma, obj.ind_p_vt, obj.ind_n_vt, obj.ind_gamma_d, obj.ind_delta_d, 'UniformOutput', false);
            else
                grab = @(l, ln, lp ,yg, g, pvt, nvt, gd,dd) vertcat(obj.w(l), obj.w(ln), obj.w(lp), obj.w(yg), obj.w(g), obj.w(pvt), obj.w(nvt), obj.w(gd), obj.w(dd));
                cross_comp_cont_0 = cellfun(grab, obj.ind_lam, obj.ind_lambda_n, obj.ind_lambda_p,...
                    obj.ind_y_gap, obj.ind_gamma, obj.ind_p_vt, obj.ind_n_vt, obj.ind_gamma_d, obj.ind_delta_d, 'UniformOutput', false);
            end

        end

        function cross_comp_cont_1 = get.cross_comp_cont_1(obj)
            import casadi.*
            if obj.model.friction_exists && obj.settings.friction_model == 'Conic' && obj.settings.conic_model_switch_handling == 'Abs'
                grab = @(pvt) vertcat(obj.w(pvt));
                cross_comp_cont_1 = cellfun(grab, obj.ind_p_vt, 'UniformOutput', false);
            else
                grab = @(pvt) [];
                cross_comp_cont_1 = cellfun(grab, obj.ind_p_vt, 'UniformOutput', false);
            end 
        end

        function cross_comp_cont_2 = get.cross_comp_cont_2(obj)
            import casadi.*
            if obj.model.friction_exists && obj.settings.friction_model == 'Conic' && obj.settings.conic_model_switch_handling == 'Abs'
                grab = @(nvt) vertcat(obj.w(nvt));
                cross_comp_cont_2 = cellfun(grab, obj.ind_n_vt, 'UniformOutput', false);
            else
                grab = @(nvt) [];
                cross_comp_cont_2 = cellfun(grab,  obj.ind_n_vt, 'UniformOutput', false);
            end
        end



        function x = get.x(obj)
            x= cellfun(@(x) obj.w(x), obj.ind_x, 'UniformOutput', false);
        end

        function y_gap = get.y_gap(obj)
            y_gap = cellfun(@(y_gap) obj.w(y_gap), obj.ind_y_gap, 'UniformOutput', false);
        end

        function gamma = get.gamma(obj)
            gamma= cellfun(@(gamma) obj.w(gamma), obj.ind_gamma, 'UniformOutput', false);
        end

        function gamma_d = get.gamma_d(obj)
            gamma_d = cellfun(@(gamma_d) obj.w(gamma_d), obj.ind_gamma_d, 'UniformOutput', false);
        end

        function delta_d = get.delta_d(obj)
            delta_d = cellfun(@(delta_d) obj.w(delta_d), obj.ind_delta_d, 'UniformOutput', false);
        end

        function p_vt = get.p_vt(obj)
            p_vt = cellfun(@(p_vt) obj.w(p_vt), obj.ind_p_vt, 'UniformOutput', false);
        end

        function n_vt = get.n_vt(obj)
            n_vt = cellfun(@(n_vt) obj.w(n_vt), obj.ind_n_vt, 'UniformOutput', false);
        end
    end
end
