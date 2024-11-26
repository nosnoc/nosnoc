classdef PDSObjects < nosnoc.dcs.Base
% This is a class that represents a Gradient complementarity system where the gap functions are
% defined implicitly via a convex optimization problem. :cite:`Tracy2023`.
%   
% It is defined by the dynamics:
%
% .. math::
%     :nowrap:
%
%     \begin{align*}
%        \dot{x} &= f(x) + \nabla c(x)\lambda \\
%        0 &\le c(x) \perp \lambda \ge 0
%     \end{align*}
%
% Where $c(x)$ is defined by:
%
% % .. math::
%     :nowrap:
%
%     \begin{align*}
%       c(x) = &\min_\alpha \alpha - 1\\
%              & \textrm{s.t}\\
%              & \quad x\in \mathcal{S}_1(\alpha)
%              & \quad x\in \mathcal{S}_2(\alpha)
%     \end{align*}
    properties
        lambda % casadi.SX|casadi.MX: Lagrange multipliers for the gradient complementarity system.
        c_lift % casadi.SX|casadi.MX: Lift variables for the gap functions in :class:`nosnoc.model.Pds`.

        x_dot_lift
        g_friction
        
        f_x % casadi.SX or casadi.MX: Right hand side of the gradient complementarity system

        dims % struct: dimensions struct TODO document members

        f_unconstrained_fun % casadi.Function: Unconstrained dynamics.
        nabla_c_fun % casadi.Function: Gradient of the gap functions.
        c_fun % casadi.Function: Gap Functions.
        g_c_lift_fun % casadi.Function: Lifting constraints for the gap functions in :class:`nosnoc.model.Pds`.

        g_kkt_fun
        g_d_fun
        g_tangent_fun
        g_friction_fun
        G_friction_fun
        H_friction_fun
        normal_fun

    end

    methods
        function obj = PDSObjects(model)
            obj.model = model;
            obj.dims = model.dims;
        end
        
        function generate_variables(obj, opts)
            import casadi.*
            model = obj.model;
            dims = obj.dims;
            
            if opts.gcs_lift_gap_functions
                dims.n_c_lift = dims.n_c;
            else
                dims.n_c_lift = 0;
            end
            obj.c_lift = define_casadi_symbolic(opts.casadi_symbolic_mode,'c_lift',dims.n_c_lift);
            obj.lambda = model.lambda;
            % All variables created at model stage for ease of use
            % TODO: Perhaps this should be moved but for now it is OK
            obj.dims = dims;
        end

        function generate_equations(obj, opts)
            import casadi.*
            model = obj.model;
            dims = obj.dims;

            dims = obj.dims;

            obj.f_x = [];

            for ii=1:length(model.objects)
                obj.x_dot_lift = vertcat(obj.x_dot_lift, model.objects(ii).x_dot_lift);
                obj.f_x = vertcat(obj.f_x, model.objects(ii).x_dot);
            end
            
            if opts.gcs_lift_gap_functions
                g_c_lift = obj.c_lift - model.c;
                c = obj.c_lift;
            else
                g_c_lift = [];
                c = model.c;
            end

            if 1 % TODO(@anton) add option for lifting x_dot
                fun_opts.allow_free = true;
                g_friction_sub_fun = Function('g_friction_sub', {obj.x_dot_lift}, {model.g_friction}, fun_opts);
                obj.g_friction = g_friction_sub_fun(obj.f_x);
            end

            dims.n_g_d = size(model.g_d, 1);
            obj.dims = dims;

            obj.g_kkt_fun = Function('g_kkt', {model.x, model.mu, model.p_d, model.y1_d, model.y2_d}, {model.g_kkt});
            obj.g_d_fun = Function('g_d', {model.x, model.alpha, model.p_d, model.y1_d, model.y2_d}, {model.g_d});
            obj.normal_fun = Function('normal', {model.x, model.mu, model.p_d, model.p}, {model.normal});
            obj.g_friction_fun = Function('g_friction', {model.x, model.u, model.lambda, model.lambda_t, model.mu, model.p_d, model.normal_lift, model.tangent, model.v_t, model.p}, {obj.g_friction});
            obj.G_friction_fun = Function('G_friction', {model.lambda, model.lambda_t, model.v_t, model.gamma_f}, {model.G_friction});
            obj.H_friction_fun = Function('H_friction', {model.lambda_t, model.gamma_f}, {model.H_friction});
            obj.f_x_fun = Function('f_x', {model.x, model.z, model.u, model.lambda, model.lambda_t, model.tangent, model.mu, model.p_d, model.v_global, model.p}, {obj.f_x});
            obj.c_fun = Function('c_fun', {model.x, obj.c_lift, model.alpha, model.z, model.v_global, model.p}, {c});
            obj.g_c_lift_fun = Function('c_fun', {model.x, obj.c_lift, model.alpha, model.z, model.v_global, model.p}, {g_c_lift});
            obj.g_z_fun = Function('g_z', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_z});
            obj.g_path_fun = Function('g_path', {model.x, model.z, model.u, model.v_global, model.p}, {model.g_path});
            obj.g_terminal_fun = Function('g_terminal', {model.x, model.z, model.v_global, model.p_global}, {model.g_terminal});
            obj.f_q_fun = Function('f_q', {model.x, model.z, model.u, model.v_global, model.p}, {model.f_q});
            obj.f_q_T_fun = Function('f_q_T', {model.x, model.z, model.v_global, model.p}, {model.f_q_T});
        end
    end

    methods(Access=protected)
        function propgrp = getPropertyGroups(obj)
            propgrp = getPropertyGroups@nosnoc.dcs.Base(obj);
            group_title = 'Variables';
            var_list = struct;
            var_list.x = obj.model.x;
            if ~any(size(obj.model.u) == 0)
                var_list.u = obj.model.u;
            end
            if ~any(size(obj.model.z) == 0)
                var_list.z = obj.model.z;
            end
            var_list.lambda = obj.lambda;
            if ~any(size(obj.c_lift) == 0)
                var_list.c_lift = obj.c_lift;
            end
            propgrp(2) = matlab.mixin.util.PropertyGroup(var_list, group_title);
        end
    end
end
