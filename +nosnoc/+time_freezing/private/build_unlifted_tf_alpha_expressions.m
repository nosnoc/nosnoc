function [alpha_ode, alpha_aux] = build_unlifted_tf_alpha_expressions(alpha_q, alpha_qv, alpha_v_normal, alpha_v_tangent, cls_model, opts, dims)
    import casadi.*
    alpha_ode = 1;
    alpha_aux = SX(zeros(dims.n_aux ,1));
    for ii = 1:dims.n_contacts
        if opts.time_freezing_nonsmooth_switching_fun
            alpha_ode = alpha_ode*alpha_qv(ii);
            if cls_model.friction_exists
                alpha_aux(ii) = (1-alpha_qv(ii))*(alpha_v_tangent(ii));
                alpha_aux(dims.n_contacts+ii) = (1-alpha_qv(ii))*(1-alpha_v_tangent(ii));
            else
                alpha_aux(ii)=(1-alpha_qv(ii));
            end
        else
            alpha_ode = alpha_ode*(alpha_q(ii)+alpha_v_normal(ii)-alpha_q(ii)*alpha_v_normal(ii));
            if cls_model.friction_exists
                alpha_aux(ii) = (1-alpha_q(ii))*(1-alpha_v_normal(ii))*(alpha_v_tangent(ii));
                alpha_aux(dims.n_contacts+ii) = (1-alpha_q(ii))*(1-alpha_v_normal(ii))*(1-alpha_v_tangent(ii));
            else
                alpha_aux(ii) = (1-alpha_q(ii))*(1-alpha_v_normal(ii));
            end
        end
    end
end
