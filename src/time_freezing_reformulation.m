function [model,settings] = time_freezing_reformulation(model,settings)
import casadi.*
%% Load settings and model details
unfold_struct(model,'caller');
time_freezing = settings.time_freezing;
time_freezing_model_exists = 0;
% TODO shift all model updates to the end of the file
% TODO add to all functions model. prefix to get rid of unfold_struct
inv_M_once = 0; % experimental option, 
    % check is the model partially filled
        if exist('F')
            fprintf('nosnoc: model.F provided, the automated model reformulation will be not performed. \n')
            time_freezing_model_exists = 1;
        else
            time_freezing_model_exists = 0;
        end

    %% Check is there a gap gunctions
    if ~exist('f_c')
        error('nosnoc: Please provide the gap functions model.f_c.')
    end
    n_contacts = length(f_c);

    % check dimensions of contacts
    if ~exist('n_dim_contact')
        error('nosnoc: Please n_dim_contact, dimension of tangent space at contact (1, 2 or 3)')
    end

    %% Model Checks
    % Check does coeffiecnt of restituion and friction exist and do they have appropiate values, check if f and M exist.
    if ~time_freezing_model_exists
        if ~exist('e')
            error('nosnoc:  Please provide a coefficient of restitution via model.e')
        end

        if abs(1-e)>1 || e<0
            error('nosnoc: the coefficient of restitution e should be in [0,1].')
        end

        if ~exist('mu')
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
        % parameter for auxiliary dynamics
        if ~exist('model.a_n') && ~exist('a_n')
            a_n  = 100;
        end
        % no friction for one dimensional models
        if n_dim_contact < 2
            % no friction for a point as there is no tangentional direction
            mu = 0;
            friction_exists  = 0;
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
        if ~exist('q','var') && ~exist('v','var')
            q = x(1:n_q);
            v = x(n_q+1:2*n_q);
        end
        


        % check model function
        if ~exist('f_v','var')
            error('nosnoc: the function f_v (collecting all generalized forces), in dv/dt =  f_v(q,v,u) + ... is not provided in model.');
        end

        if ~exist('f_gravity','var')
            if ~exist('gravity','var')
                f_gravity = zeros(n_q,1);
            else
                f_gravity = gravity;
            end
        end

        % Check intertia matrix
        if ~exist('M','var')
            fprintf('nosnoc: Inertia matrix M is not provided, set to default value: M = eye(n_q). \n')
            M = eye(n_q);
            invM = inv(M);
            model.M = M;
        else
            invM = inv(M);
        end

        n_quad  = 0;
        if settings.time_freezing_quadrature_state
            % define quadrature state
            L = define_casadi_symbolic(casadi_symbolic_mode,'L',1);
            % update lower and upper bounds of lbx and ubx
            if exist('lbx')
                model.lbx = [model.lbx;-inf];
            end
            if exist('ubx')
                model.ubx = [model.ubx;inf];
            end
            x = [x;L];
            model.x = x;
            model.x0 = [model.x0;0];
            f = [f;f_q];
            model.f = f;
            model.f_q = 0;
            if exist('f_q_T','var')
                model.f_q_T  = model.f_q_T + L;
            else
                model.f_q_T  = L;
            end
            n_quad = 1;
        end
        model.n_quad = n_quad;
        %%  Normal Contact Jacobian
        try
            J_normal = model.J_normal;
            J_normal_exists = 1;
        catch
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
        %% Tangent Contact Jacobian
        try
            J_tangent = model.J_tangent;
            J_tangent_exists = 1;
        catch
            J_tangent_exists = 0;
        end

        if J_tangent_exists
            if size(J_tangent,1)~=n_q && size(J_tangent,2)~=n_contacts*(n_dim_contact-1)
                fprintf('nosnoc: J_tangentshould be %d x %d matrix.\n',n_q,n_contacts*(n_dim_contact-1));
                error('J_tangent has the wrong size.')
            end
            J_tangent_exists = 1;
        else
            fprintf('nosnoc: tangent Jacobian not provided.\n');
            J_tangent_exists = 0;
        end

        % Check is tangent jacobian avilable
        if any(mu>0) && ~J_tangent_exists
            error('Please provide the tangential contact jacobian J_tangent (the matrix contining all vectors that span the tangent spaces at contact points) or set mu = 0.');
        end

        %% Clock state and dimensions
        if ~mod(n_x,2)
            % uneven number of states = it is assumed that the clock state
            % is defined.
            t = define_casadi_symbolic(casadi_symbolic_mode,'t',1);
            % update lower and upper bounds of lbx and ubx
            if exist('lbx')
                model.lbx = [model.lbx;-inf];
            end
            if exist('ubx')
                model.ubx = [model.ubx;inf];
            end
            x = [x;t];
            model.x = x;
            model.x0 = [model.x0;0];
        end

        % compute normal vector
        if ~exist('nabla_q_f_c','var')
            nabla_q_f_c = f_c.jacobian(q)';
        end

        if is_zero(nabla_q_f_c)
            error('nosnoc: The normal vector should have at least one non-zero entry.')
        end
        F = {}; c = {};  S = {};

        %% Time-freezing reformulation
        if e == 0
            settings.time_freezing_inelastic = 1; % flag tha inealstic time-freezing is using (for hand crafted lifting)
            settings.pss_mode = 'Step'; % currently only with step mode supported
            %% unconstrained dynamcis with clock state
            inv_M = inv(M);
             f_ode = [v;...
                    inv_M*f_v;
                     1];
            %% Auxiliary dynamics
                % where to use invM, in every aux dyn or only at the end
                if inv_M_once
                    inv_M_aux = eye(n_q);
                else
                    inv_M_aux = inv_M;
                end
            % create auxiliary dynamics
            switch n_contacts
                case 1
                    f_aux_n1 = [zeros(n_q,1);invM*nabla_q_f_c*a_n;zeros(n_quad+1,1)];
                case 2
                    if n_dim_contact >2
                        error('Multiple contacts currently only for planar models avilable.')
                    end
                    f_aux_n1 = [zeros(n_q,1);invM*nabla_q_f_c(:,1)*a_n;zeros(n_quad+1,1)];
                    f_aux_n2 = [zeros(n_q,1);invM*nabla_q_f_c(:,2)*a_n;zeros(n_quad+1,1)];
                otherwise
                    error('Currently only up to two unilateral constraint in the planar case suported.')
            end

            % compute tangents and frictional auxiliary dynamics
            if friction_exists
               
                % create auxiliary dynamics
                a_t = mu*a_n;
                switch n_dim_contact
                    case 2              
                        f_aux_t1 = [zeros(n_q,1);invM*J_tangent(:,1)*a_t;zeros(n_quad+1,1)];
                        f_aux_t2 = -f_aux_t1;
                    case 3
                        eps_tangential = 1e-8;
                        T_q = [tangent1 tangent2]; % spans tangent space
                        v_tan = T_q'*v;
                        v_tan_norm = norm(v_tan+1e-16); % avoid dividing by zero at initalisation.
                        % auxiliary dynamics pushing towards origin
                        f_aux_t1 = [zeros(n_q,1);-invM*T_q*(a_t*v_tan/(v_tan_norm));zeros(n_quad+1,1)];
                        % relaxed pushing toward relaxed circle around origin
                        f_aux_t2 = [zeros(n_q,1);invM*T_q*v_tan;zeros(n_quad+1,1)];
                end
            end
            % switching functions
            c1 = f_c;
            c2 = nabla_q_f_c'*v;
            if friction_exists
                switch n_dim_contact
                    case 2
                        c3 = J_tangent(:,1)'*v;
                    case 3
                        c3 = eps_tangential-v_tan_norm;
                end
            else
                c3 = [];
            end

            if ~friction_exists
                switch n_contacts
                    case 1
                        F{1} = [f_ode f_aux_n1];
                        S{1} = [1 0;-1 1;-1 -1];
                    case 2
                        F{1} = [f_ode f_aux_n1 f_aux_n2 f_aux_n1+f_aux_n2];
                        S{1} = eye(4);
                    otherwise
                        error('not implemented.')
                end
            else
                switch n_contacts
                    case 1
                        F{1} = [f_ode,f_aux_n1+f_aux_t1,f_aux_n1+f_aux_t2];
                        S{1} = [1 0 0;-1 -1 -1;-1 -1 1];
                    case 2
                        F{1} = [f_ode f_aux_n1 f_aux_n2 f_aux_n1+f_aux_n];
                        S{1} = eye(4);
                        F{1} = [f_ode,f_aux_n1+f_aux_t1,f_aux_n1+f_aux_t2,f_aux_n2+f_aux_t3,f_aux_n2+f_aux_t4, f_aux_n1+f_aux_n2];
                end
            end

            c{1} = [c1;c2;c3];
            % store results and flag that model is created
            model.F = F; 
            model.c = c;
            model.S = S;
            time_freezing_model_exists = 1;
        else
            % elastic
            if ~exist('k_aux')
                k_aux = 10;
                if settings.print_level > 1
                    fprintf('nosnoc: Setting default value for k_aux = 10.\n')
                end
            end
            temp1 =2*abs(log(e));
            temp2 = k_aux/(pi^2+log(e)^2);
            c_aux = temp1/sqrt(temp2);
            %                 c_aux = 0.211989;
            %             f_aux_n = [0;v(2);0;-k_aux*q(2)-c_aux*v(2);0];
            K = [0 1;-k_aux -c_aux];
            N  = [nabla_q_f_c zeros(n_q,1);...
                zeros(n_q,1) invM*nabla_q_f_c];
            f_aux_n1 = N*K*N'*[q;v];
            f_aux_n1 = [f_aux_n1;zeros(n_quad+1,1)];
            f = [v;f;1];
            % updated with clock state
            model.f = f;
            model.F = [f, f_aux_n1];
            model.S = [1; -1];
            time_freezing_model_exists = 1;
        end
    end

% fprintf('nosnoc: No action was done. Consider setting settings.time_freezing = 1, if calling this function.\n')

%% Settings updates
settings.time_freezing_model_exists = time_freezing_model_exists;
settings.friction_exists  = friction_exists;

%% Model updates
model.n_q = n_q;
model.q = q;
model.v = v;
model.n_contacts = n_contacts;
end