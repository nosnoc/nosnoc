function [alpha, alpha_q, alpha_qv, alpha_v_normal, alpha_v_tangent] = generate_tf_inelastic_alphas(cls_model, tf_model, opts, dims)
    alpha = define_casadi_symbolic(opts.casadi_symbolic_mode,'alpha', length(tf_model.c));
    alpha_q = [];
    alpha_qv = [];
    alpha_v_normal = [];
    alpha_v_tangent = [];
    if ~opts.time_freezing_nonsmooth_switching_fun
        alpha_q = alpha(1:dims.n_contacts);
        alpha_v_normal = alpha(dims.n_contacts+1:2*dims.n_contacts);
        if cls_model.friction_exists
            alpha_v_tangent = alpha(2*dims.n_contacts+1:end);
        end
    else
        alpha_qv = alpha(1:dims.n_contacts);
        if cls_model.friction_exists
            alpha_v_tangent = alpha(dims.n_contacts+1:end);
        end
    end
end
