function pss_model = cls_inelastic(cls_model, opts)
    import casadi.*
    dims = cls_model.dims;
    pss_model = nosnoc.model.Pss();
    pss_model.dims = dims;

    % transfer all common properties to pss_model
    common_props = properties('nosnoc.model.Base');
    for ii=1:length(common_props)
        prop = common_props{ii};
        pss_model.(prop) = cls_model.(prop);
    end
    
    dims.n_contacts = length(cls_model.f_c);
    % check dimensions of contacts
    if isempty(dims.n_dim_contact)
        warning('nosnoc: Please set n_dim_contact, dimension of tangent space at contact (1, 2)')
        dims.n_dim_contact = 1;
    end
    
    % quadrature state
    dims.n_quad = 0;
    f_quad = [];
    if opts.time_freezing_quadrature_state
        % define quadrature state
        L = define_casadi_symbolic(opts.casadi_symbolic_mode,'L',1);
        pss_model.lbx = [pss_model.lbx;-inf];
        pss_model.ubx = [pss_model.ubx;inf];
        pss_model.x = [pss_model.x;L];
        pss_model.x0 = [pss_model.x0;0];
        f_quad = pss_model.f_q;
        pss_model.f_q = 0;
        if ~isempty(pss_model.f_q_T)
            pss_model.f_q_T  = pss_model.f_q_T + L;
        else
            pss_model.f_q_T = L;
        end
        dims.n_quad = 1;
    end
    
    % uneven number of states = it is assumed that the clock state is defined.
    t = define_casadi_symbolic(opts.casadi_symbolic_mode,'t',1);
    % update lower and upper bounds of lbx and ubx
    pss_model.lbx = [pss_model.lbx;-inf];
    pss_model.ubx = [pss_model.ubx;inf];
    pss_model.x = [pss_model.x;t];
    pss_model.x0 = [pss_model.x0;0];
    
    % normal and tangential velocities
    eps_t = opts.eps_t;
    v_normal = cls_model.J_normal'*cls_model.v;
    v_tangent = [];
    v_tangent_norms = [];
    if cls_model.friction_exists
        if dims.n_dim_contact == 1
            v_tangent = (cls_model.J_tangent'*cls_model.v)';
        else
            v_tangent = cls_model.J_tangent'*cls_model.v;
            v_tangent = reshape(v_tangent,2,dims.n_contacts); % 2 x n_c , the columns are the tangential velocities of the contact points
        end
        v_tangent_norms = [];
        for ii = 1:dims.n_contacts
            v_tangent_norms = [v_tangent_norms;norm(v_tangent(:,ii))];
        end
    else
        v_tangent  = [];
    end
    
    % parameter for auxiliary dynamics
    a_n = opts.a_n;
    %% Time-freezing reformulation
    % unconstrained dynamics with clock state
    f_ode = [cls_model.v;...
        cls_model.invM*cls_model.f_v;
        f_quad;
        1];

    % Auxiliary dynamics
    inv_M_aux = cls_model.invM;
    inv_M_ext = eye(dims.n_x+1+dims.n_quad);
    f_aux_pos = []; % matrix with all auxiliary tangent dynamics
    f_aux_neg = [];
    % time freezing dynamics
    if opts.stabilizing_q_dynamics
        f_q_dynamics = -opts.kappa_stabilizing_q_dynamics*cls_model.J_normal*diag(cls_model.f_c);
    else
        f_q_dynamics = zeros(dims.n_q,dims.n_contacts);
    end
    f_aux_normal = [f_q_dynamics;inv_M_aux*cls_model.J_normal*a_n;zeros(1+dims.n_quad, dims.n_contacts)];

    if opts.time_freezing_nonsmooth_switching_fun
        pss_model.c = [max_smooth_fun(cls_model.f_c,v_normal,0);v_tangent];    
    else
        if dims.n_dim_contact == 1
            pss_model.c = [cls_model.f_c;v_normal;v_tangent'];
        else
            pss_model.c = [cls_model.f_c;v_normal;v_tangent_norms-eps_t];
        end
    end

    for ii = 1:dims.n_contacts
        if cls_model.friction_exists && cls_model.mu(ii)>0
            % auxiliary tangent;
            if dims.n_dim_contact == 1
                v_tangent_ii = cls_model.J_tangent(:,ii)'*cls_model.v;
                f_aux_pos_ii = [f_q_dynamics(:,ii) ;inv_M_aux*(cls_model.J_normal(:,ii)-cls_model.J_tangent(:,ii)*(cls_model.mu(ii)))*a_n;zeros(1+dims.n_quad,1)]; % for v>0
                f_aux_neg_ii = [f_q_dynamics(:,ii) ;inv_M_aux*(cls_model.J_normal(:,ii)+cls_model.J_tangent(:,ii)*(cls_model.mu(ii)))*a_n;zeros(1+dims.n_quad,1)]; % for v<0
            else
                v_tangent_ii = v_tangent(:,ii);
                f_aux_pos_ii = [f_q_dynamics(:,ii);inv_M_aux*(cls_model.J_normal(:,ii)*a_n-cls_model.J_tangent(:,ii*2-1:ii*2)*cls_model.mu(ii)*a_n*v_tangent_ii/norm(v_tangent_ii+eps_t));zeros(1+dims.n_quad,1)]; % for v>0
                f_aux_neg_ii = [f_q_dynamics(:,ii);inv_M_aux*(cls_model.J_normal(:,ii)*a_n+cls_model.J_tangent(:,ii*2-1:ii*2)*cls_model.mu(ii)*a_n*v_tangent_ii/norm(v_tangent_ii+eps_t));zeros(1+dims.n_quad,1)]; % for v>0
            end
            f_aux_pos = [f_aux_pos,f_aux_pos_ii];
            f_aux_neg = [f_aux_neg,f_aux_neg_ii];
        end
    end
    if cls_model.friction_exists
        f_aux = [f_aux_pos,f_aux_neg];
        dims.n_aux = 2*dims.n_contacts;
    else
        f_aux = f_aux_normal;
        dims.n_aux = dims.n_contacts;
    end

    % Build logical functions because its easier than hand picking logic.
    [alpha, alpha_q, alpha_qv, alpha_v_normal, alpha_v_tangent] = generate_tf_inelastic_alphas(cls_model, pss_model, opts, dims);
    [alpha_ode, alpha_aux] = build_unlifted_tf_alpha_expressions(alpha_q, alpha_qv, alpha_v_normal, alpha_v_tangent, cls_model, opts, dims);
    % This function when passed 0,1 arguments returns 1 for which f should be used for the
    % region defined by the arguments.
    indicator = Function('alpha_indicator', {alpha}, {[alpha_ode; alpha_aux]});

    % This generates the dense S matrix.
    S = dec2bin(0:2^((2+cls_model.friction_exists)*dims.n_contacts)-1) - '0';
    S = flip(S,1);
    indicators = full(indicator(S')');

    % Generate vector of all f.
    f_all = [f_ode (inv_M_ext*f_aux)];
    F = {};

    % Populate F by going through each row of S and selecting the correct f.
    for ii=1:size(S,1)
        ind = find(indicators(ii,:));
        f_i = f_all(:,ind);
        F{ii} = f_i;
    end
    
    pss_model.F = horzcat(F{:});

    % Correct S to have -1 where the step function is zero for the corresponding c(x).
    S(~S) = -1;
    pss_model.S = S;
    
    pss_model.dims = dims;
end
