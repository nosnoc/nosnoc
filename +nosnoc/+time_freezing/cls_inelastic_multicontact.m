function heaviside_model = cls_inelastic_multicontact(cls_model, opts)
% Implements the custom lifting formulation that takes advantage of the inelastic time freezing formulation.
% This allows for a non-exponential number of regions as a function of the number of contacts.
% This is done by taking advantage of the shared free-flight region and the joint auxiliary dynamics:
% $$F_{\mathrm{aux}} = \sum_{i=1}^{n_c} F_{\mathrm{aux}}^i(x)$$
    import casadi.*
    dims = cls_model.dims;
    heaviside_model = nosnoc.model.Heaviside();
    heaviside_model.dims = dims;

    % transfer all common properties to heaviside_model
    common_props = properties('nosnoc.model.Base');
    for ii=1:length(common_props)
        prop = common_props{ii};
        heaviside_model.(prop) = cls_model.(prop);
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
        heaviside_model.lbx = [heaviside_model.lbx;-inf];
        heaviside_model.ubx = [heaviside_model.ubx;inf];
        heaviside_model.x = [heaviside_model.x;L];
        heaviside_model.x0 = [heaviside_model.x0;0];
        f_quad = heaviside_model.f_q;
        heaviside_model.f_q = 0;
        if ~isempty(heaviside_model.f_q_T)
            heaviside_model.f_q_T  = heaviside_model.f_q_T + L;
        else
            heaviside_model.f_q_T = L;
        end
        dims.n_quad = 1;
    end
    
    % uneven number of states = it is assumed that the clock state is defined.
    t = define_casadi_symbolic(opts.casadi_symbolic_mode,'t',1);
    % update lower and upper bounds of lbx and ubx
    heaviside_model.lbx = [heaviside_model.lbx;-inf];
    heaviside_model.ubx = [heaviside_model.ubx;inf];
    heaviside_model.x = [heaviside_model.x;t];
    heaviside_model.x0 = [heaviside_model.x0;0];
    
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
    f_aux_normal = [f_q_dynamics;inv_M_aux*cls_model.J_normal*a_n;zeros(1+dims.n_quad,dims.n_contacts)];

    if opts.time_freezing_nonsmooth_switching_fun
        sigma = SX.sym('sigma',1);
        a = SX.sym('a',1);
        b = SX.sym('b',1);
        f_natural_residual = 0.5*(b+a+sqrt((b-a+sigma)^2));
        % f_natural_residual = max(a,b);
        max_smooth_fun = Function('max_smooth_fun',{a,b,sigma},{f_natural_residual});
        heaviside_model.c = [max_smooth_fun(cls_model.f_c,v_normal,0);v_tangent];    
    else
        if dims.n_dim_contact == 1
            heaviside_model.c = [cls_model.f_c;v_normal;v_tangent'];
        else
            heaviside_model.c = [cls_model.f_c;v_normal;v_tangent_norms-eps_t];
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

    [alpha, alpha_q, alpha_qv, alpha_v_normal, alpha_v_tangent] = generate_tf_inelastic_alphas(cls_model, heaviside_model, opts, dims);
    theta = define_casadi_symbolic(opts.casadi_symbolic_mode,'theta',dims.n_aux+1);

    % lifting variables
    % empty expressions for initalization
    beta_bilinear_ode = []; % for lifting bilinear terms in free flight dynamics multiplire
    beta_bilinear_aux = []; % for lifiting bilinear terms appearing in aux. dynamics mutiplieres
    beta_prod = []; % for lifting the multi affine term defineng the overall free flight dynamics multpliers
                    % expressions for lifting
    beta_bilinear_ode_expr = [];
    beta_bilinear_aux_expr = [];
    beta_prod_expr = [];
    beta_prod_expr_guess = []; % extra expresion to make depend only on alpha (the one above depens on both and alpha and beta) - needed for eval. of inital guess

    alpha_ode = 1;
    if opts.casadi_symbolic_mode == "casadi.MX"
        alpha_aux = MX(zeros(dims.n_aux, 1));
    else
        alpha_aux = SX(zeros(dims.n_aux, 1));
    end
    % lift bilinear terms in product terms for free flight ode % (alpha_q*alpha_v)
    if ~opts.time_freezing_nonsmooth_switching_fun
        beta_bilinear_ode = define_casadi_symbolic(opts.casadi_symbolic_mode,'beta_bilinear_ode',dims.n_contacts);
        beta_bilinear_ode_expr = eval([opts.casadi_symbolic_mode '.zeros(' num2str(dims.n_contacts) ',1);']);
        if cls_model.friction_exists
            % lift bilinear terms defining aux dynamics (1-alpha_q)*(1-alpha_v)
            beta_bilinear_aux = define_casadi_symbolic(opts.casadi_symbolic_mode,'beta_bilinear_aux',dims.n_contacts);
            beta_bilinear_aux_expr = eval([opts.casadi_symbolic_mode '.zeros(' num2str(dims.n_contacts) ',1);']);
        end
    end
    if dims.n_contacts > 2
        beta_prod = define_casadi_symbolic(opts.casadi_symbolic_mode,'beta_bilinear',dims.n_contacts-2);
        beta_prod_expr = eval([opts.casadi_symbolic_mode '.zeros(' num2str(dims.n_contacts) '-2,1);']);
        beta_prod_expr_guess = eval([opts.casadi_symbolic_mode '.zeros(' num2str(dims.n_contacts) '-2,1);']);
    end
    
    % lift and have bilinear terms
    if opts.time_freezing_nonsmooth_switching_fun
        if dims.n_contacts <= 2
            for ii = 1:dims.n_contacts
                alpha_ode = alpha_ode*alpha_qv(ii);
            end
        else
            beta_prod_expr(1) = (alpha_qv(1))*(alpha_qv(2));
            beta_prod_expr_guess(1) = (alpha_qv(1))*(alpha_qv(2));
            % lifting terms in between
            for ii = 3:dims.n_contacts-1
                beta_prod_expr(ii-1) = beta_prod(ii-2)*(alpha_qv(ii)); % beta_{i} = beta{i-1}*(prod_term_i+1}
                beta_prod_expr_guess(ii-1) = beta_prod_expr_guess(ii-2)*(alpha_qv(ii)); % this is to have an expression depending only on alpha for the inital guess eval
            end
            alpha_ode = beta_prod(end)*(alpha_qv(dims.n_contacts)); % last lifting term;
        end
        % lifting of aux dyn multiplier expressions
        for ii = 1:dims.n_contacts
            if cls_model.friction_exists
                alpha_aux(ii) = (1-alpha_qv(ii))*(alpha_v_tangent(ii));
                alpha_aux(dims.n_contacts+ii) = (1-alpha_qv(ii))*(1-alpha_v_tangent(ii));
            else
                alpha_aux(ii)=(1-alpha_qv(ii));
            end
        end
    else
        % with two smooth switching functions
        beta_bilinear_ode_expr = alpha_q.*alpha_v_normal;
        if cls_model.friction_exists
            beta_bilinear_aux_expr = (1-alpha_q).*(1-alpha_v_normal);
        end

        if dims.n_contacts <= 2
            % here no lifting of product terms
            for ii = 1:dims.n_contacts
                alpha_ode = alpha_ode*(alpha_q(ii)+alpha_v_normal(ii)-beta_bilinear_ode(ii));
            end
        else
            % here lifting of product terms
            beta_prod_expr(1) = [(alpha_q(1)+alpha_v_normal(1)-beta_bilinear_ode(1))*(alpha_q(2)+alpha_v_normal(2)-beta_bilinear_ode(2))]; % first lifting terms
                                                                                                                                                          % lifting terms in between
            for ii = 3:dims.n_contacts-1
                beta_prod_expr(ii-1) = beta_prod(ii-2)*(alpha_q(ii)+alpha_v_normal(ii)-beta_bilinear_ode(ii));
                beta_prod_expr_guess(ii-1) = beta_prod_expr_guess(ii-2)*(alpha_q(ii)+alpha_v_normal(ii)-beta_bilinear_ode(ii));
            end
            % last lifting term;
            alpha_ode = beta_prod(end)*(alpha_q(dims.n_contacts)+alpha_v_normal(dims.n_contacts)-beta_bilinear_ode(dims.n_contacts));
        end
        % lifting of aux dyn multiplier expressions
        for ii = 1:dims.n_contacts
            if cls_model.friction_exists
                alpha_aux(ii) = beta_bilinear_aux(ii)*(alpha_v_tangent(ii));
                alpha_aux(dims.n_contacts+ii) = beta_bilinear_aux(ii)*(1-alpha_v_tangent(ii));
            else
                alpha_aux(ii) = (1-alpha_q(ii))*(1-alpha_v_normal(ii));
            end
        end
    end

    beta = [beta_bilinear_ode;
        beta_bilinear_aux;
        beta_prod];

    try
        beta0_fun = Function('beta0', {alpha}, {[beta_bilinear_ode_expr; beta_bilinear_aux_expr; beta_prod_expr]});
        theta0_fun = Function('theta0', {alpha,beta}, {[alpha_ode;alpha_aux]});
        
        beta0 = full(beta0_fun(opts.initial_alpha*ones(size(alpha,1),1)));
        theta0 = full(theta0_fun(opts.initial_alpha*ones(size(alpha,1),1),beta0));
    catch e
        % TODO(@anton) Mild hack, improve later.
        w_init = [beta;theta];
        p_init = [alpha];
        g_init = [beta - [beta_bilinear_ode_expr; beta_bilinear_aux_expr; beta_prod_expr]; theta - [alpha_ode;alpha_aux]];
        nlp_init = struct('x', w_init, 'f', 0, 'g', g_init, 'p', p_init);

        opts_init.ipopt.sb = 'no';
        opts_init.ipopt.print_level = 0;
        init = nlpsol('init', 'ipopt', nlp_init, opts_init);
        
        init_res = init('p', opts.initial_alpha*ones(size(alpha,1),1), 'lbg', zeros(size(g_init)), 'ubg', zeros(size(g_init)));
        
        beta0 = full(init_res.x(1:size(beta, 1)));
        theta0 = full(init_res.x((size(beta, 1)+1):end));
    end
    heaviside_model.z = [theta;beta];
    heaviside_model.z0 = [theta0;beta0];
    heaviside_model.g_z = [theta-[alpha_ode;alpha_aux];
        beta - [beta_bilinear_ode_expr; beta_bilinear_aux_expr; beta_prod_expr]];
    
    f_all = [f_ode (inv_M_ext*f_aux)];
    heaviside_model.f_x = f_all*theta;
    heaviside_model.alpha = alpha;

    pss_model.dims = dims;
end
