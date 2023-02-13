function [model,settings] = time_freezing_reformulation(model,settings)
import casadi.*
%% Load settings and model details
unfold_struct(model,'caller');
time_freezing = settings.time_freezing;
time_freezing_model_exists = 0;
%% TODOS
% %shift all model updates to the end of the file
% add to all functions model. prefix to get rid of unfold_struct
% refactor time_freezing_quadrature_state = true
% check is the model partially filled
if isfield(model,'F')
    fprintf('nosnoc: model.F provided, the automated model reformulation will be not performed. \n')
    time_freezing_model_exists = 1;
else
    time_freezing_model_exists = 0;
end

%% Experimental options
inv_M_once = 0;
%% Auxiliary functions
sigma = SX.sym('sigma',1);
a = SX.sym('a',1);
b = SX.sym('b',1);
f_natural_residual = 0.5*(b+a+sqrt((b-a+sigma)^2));
% f_natural_residual = max(a,b);
max_smooth_fun = Function('max_smooth_fun',{a,b,sigma},{f_natural_residual});
%% Chek is the provided user data valid and complete
% Check is there a gap gunctions
if ~isfield(model,'f_c')
    error('nosnoc: Please provide the gap functions model.f_c.')
end
n_contacts = length(f_c);
% check dimensions of contacts
if ~isfield(model,'n_dim_contact')
    error('nosnoc: Please n_dim_contact, dimension of tangent space at contact (1, 2 or 3)')
end

if ~time_freezing_model_exists
    % coefficent of restiution
    if ~isfield(model,'e')
        error('nosnoc:  Please provide a coefficient of restitution via model.e')
    end

    if abs(1-e)>1 || e<0
        error('nosnoc: the coefficient of restitution e should be in [0,1].')
    end

    % coefficient of friction
    if ~isfield(model,'mu')
        mu = 0;
        fprintf('nosnoc: Coefficients of friction mu not provided, setting it to zero for all contacts. \n')
    end

    if length(mu(:)) ~= 1 && length(mu(:)) ~= n_contacts
        errro('nosnoc: Vector mu has to have length 1 or n_c.')
    end

    if any(mu<0)
        error('nosnoc: The coefficients of friction mu should be nonnegative.')
    end

    if any(mu)>0
        friction_exists = 1;
        if length(mu(:)) == 1
            mu = ones(n_contacts,1)*mu;
        end
    else
        friction_exists = 0;
    end

    % dimensions and state space split
    casadi_symbolic_mode = model.x(1).type_name();
    if mod(size(x,1),2)
        n_x = size(x,1);
        n_q = (n_x-1)/2;
    else
        n_x = size(x,1);
        n_q = n_x/2;
    end

    if ~isfield(model,'q') && ~isfield(model,'v')
        q = x(1:n_q);
        v = x(n_q+1:2*n_q);
    end

    % check model function
    if ~isfield(model,'f_v')
        error('nosnoc: the function f_v (collecting all generalized forces), in dv/dt =  f_v(q,v,u) + ... is not provided in model.');
    end

    % Check intertia matrix
    if ~isfield(model,'M')
        fprintf('nosnoc: Inertia matrix M is not provided, set to default value: M = eye(n_q). \n')
        M = eye(n_q);
        invM = inv(M);
    else
        invM = inv(M);
    end


    %  Normal Contact Jacobian
    if isfield(model,'J_normal')
        J_normal = model.J_normal;
        J_normal_exists = 1;
    else
        J_normal_exists = 0;
    end

    if J_normal_exists
        if size(J_normal,1)~=n_q && size(J_normal,2)~=n_contacts
            fprintf('nosnoc: J_normal should be %d x %d matrix.\n',n_q,n_contacts);
            error('nosnoc: J_normal has the wrong size.')
        end
        J_normal_exists = 1;
    else
        J_normal = f_c.jacobian(q)';
        fprintf('nosnoc: normal contact Jacobian not provided, but it is computed from the gap functions.\n');
        J_normal_exists = 1;
    end

    if is_zero(J_normal)
        error('nosnoc: The normal vector should have at least one non-zero entry.')
    end

    % Tangent Contact Jacobian
    if friction_exists
        if isfield(model,'J_tangent')
            J_tangent = model.J_tangent;
            J_tangent_exists = 1;
        else
            J_tangent_exists = 0;
        end

        if J_tangent_exists
            if size(J_tangent,1)~=n_q && size(J_tangent,2)~=n_contacts*(n_dim_contact-1)
                fprintf('nosnoc: J_tangent should be %d x %d matrix.\n',n_q,n_contacts*(n_dim_contact-1));
                error('nosnoc: J_tangent has the wrong size.')
            end
            J_tangent_exists = 1;
        else
            error('nosnoc: tangent Jacobian model.J_tangent not provided.\n');
        end
    else
        J_tangent_exists = 0;
    end

    % qudrature state
    n_quad  = 0;
    if settings.time_freezing_quadrature_state
        % define quadrature state
        L = define_casadi_symbolic(casadi_symbolic_mode,'L',1);
        if isfield(model,'lbx')
            model.lbx = [model.lbx;-inf];
        end
        if isfield(model,'ubx')
            model.ubx = [model.ubx;inf];
        end
        x = [x;L];
        model.x = x;
        model.x0 = [model.x0;0];
        f = [f;f_q];
        model.f = f;
        model.f_q = 0;
        if isfield(model,'f_q_T')
            model.f_q_T  = model.f_q_T + L;
        else
            model.f_q_T  = L;
        end
        n_quad = 1;
    end
    % Clock state and dimensions
    if ~mod(n_x,2)
        % uneven number of states = it is assumed that the clock state is defined.
        t = define_casadi_symbolic(casadi_symbolic_mode,'t',1);
        % update lower and upper bounds of lbx and ubx
        if isfield(model,'lbx')
            model.lbx = [model.lbx;-inf];
        end
        if isfield(model,'ubx')
            model.ubx = [model.ubx;inf];
        end
        x = [x;t];
        x0 = [x0;0];
    end

    % normal and tangential velocities
    eps_t = 1e-7;
    v_normal = J_normal'*v;
    if friction_exists
        if n_dim_contact == 2
            v_tangent = J_tangent'*v;
        else
            v_tangent = J_tangent'*v;
            v_tangent = reshape(v_tangent,2,n_contacts); % 2 x n_c , the columns are the tangential velocities of the contact points
           
        end
         v_tangent_norms = [];
            for ii = 1:n_contacts
                v_tangent_norms = [v_tangent_norms;norm(v_tangent(:,ii))];
            end
    else
        v_tangent  = [];
    end

  

    % parameter for auxiliary dynamics
    if ~isfield(model,'a_n')
        a_n  = 100;
    end
    %% Time-freezing reformulation
    if e == 0
        % Basic settings
        settings.time_freezing_inelastic = 1; % flag tha inealstic time-freezing is using (for hand crafted lifting)
        settings.pss_mode = 'Step'; % time freezing inelastic works better step (very inefficient with stewart)
          %% switching function
        if settings.nonsmooth_switching_fun  
            c = [max_smooth_fun(f_c,v_normal,0);v_tangent];
        %         c = max_smooth_fun(f_c,v_normal,sigma0);
        else
            c = [f_c;v_normal;v_tangent];        
        end
        %% unconstrained dynamcis with clock state
        inv_M = inv(M);
        f_ode = [v;...
            inv_M*f_v;
            1];
        
        %% Auxiliary dynamics
        % where to use invM, in every aux dyn or only at the end
        if inv_M_once
            inv_M_aux = eye(n_q);
            inv_M_ext = blkdiag(zeros(n_q),inv_M,0);
        else
            inv_M_aux = inv_M;
            inv_M_ext = eye(n_x+1);
        end
        f_aux_pos = []; % matrix wit all aux tan dyn
        f_aux_neg = [];
        % time freezing dynamics
        f_aux_normal = [zeros(n_q,n_contacts);inv_M_aux*J_normal*a_n;zeros(1,n_contacts)];
        for ii = 1:n_contacts
            if friction_exists && mu(ii)>0
                % auxiliary tangent;
                if n_dim_contact == 2
                    v_tangent_ii = J_tangent(:,ii)'*v;
                    f_aux_pos_ii = [zeros(n_q,1);inv_M_aux*(J_normal(:,ii)-J_tangent(:,ii)*(mu(ii)))*a_n;0]; % for v>0
                    f_aux_neg_ii = [zeros(n_q,1);inv_M_aux*(J_normal(:,ii)+J_tangent(:,ii)*(mu(ii)))*a_n;0]; % for v<0
                else
                    v_tangent_ii = v_tangent(:,ii);
                    f_aux_pos_ii = [zeros(n_q,1);inv_M_aux*(J_normal(:,ii)*a_n-J_tangent(:,ii*2-1:ii*2)*mu(ii)*a_n*v_tangent_ii/norm(v_tangent_ii+1e-12));0]; % for v>0
                    f_aux_neg_ii = [zeros(n_q,1);inv_M_aux*(J_normal(:,ii)*a_n+J_tangent(:,ii*2-1:ii*2)*mu(ii)*a_n*v_tangent_ii/norm(v_tangent_ii+1e-12));0]; % for v>0
                end
                f_aux_pos = [f_aux_pos,f_aux_pos_ii];
                f_aux_neg= [f_aux_neg,f_aux_neg_ii];
            end
        end
        % f_aux_normal = inv_M_aux*J_normal*a_n;
        % f_aux_tangent = inv_M_aux*J_tangent*mu(ii)*a_n;
        if friction_exists
            f_aux = [f_aux_pos,f_aux_neg];
        else
            f_aux = f_aux_normal;
        end
        F = [f_ode (inv_M_ext*f_aux)];
        S = ones(size(F,2),length(c)); % dummy value to pass error checks
        % number of auxiliary dynamicsm modes
        if friction_exists
            n_aux = 2*n_contacts;
        else
            n_aux = n_contacts;
        end
    else
        % elastic
        pss_mode = 'Step';
        if ~isfield(model,'k_aux')
            k_aux = 10;
            if settings.print_level > 1
                fprintf('nosnoc: Setting default value for k_aux = 10.\n')
            end
        end
        temp1 = 2*abs(log(e));
        temp2 = k_aux/(pi^2+log(e)^2);
        c_aux = temp1/sqrt(temp2);
        K = [0 1;-k_aux -c_aux];
        N  = [J_normal zeros(n_q,1);...
            zeros(n_q,1) invM*J_normal];
        f_aux_n1 = N*K*N'*[q;v];
        f_aux_n1 = [f_aux_n1;zeros(n_quad+1,1)];
        f_ode = [v;invM*f_v;1];
        % updated with clock state
        F = [f_ode, f_aux_n1];
        S = [1; -1];
        n_aux = 1;
        c = f_c;
    end
        time_freezing_model_exists = 1; % mark that model was created
end

%% Settings updates
settings.time_freezing_model_exists = time_freezing_model_exists;
settings.friction_exists  = friction_exists;

%% Model updates
model.n_quad = n_quad;
model.n_q = n_q;
model.n_aux = n_aux;
model.q = q;
model.v = v;
model.x = x;
model.x0 = x0;
model.M = M;
model.n_contacts = n_contacts;
model.mu = mu;
model.J_normal = J_normal;
model.F = F;
model.c = c;
model.S = S;
end