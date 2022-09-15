function [model,settings] = time_freezing_reformulation(model,settings)
import casadi.*
%% Load settings and model details
unfold_struct(model,'caller');
time_freezing = settings.time_freezing;
time_freezing_model_exists = 0;
if time_freezing
    % sanity heck
    % check is the model partially filled
    if exist('S')
        if exist('F')
            fprintf('Info on Time-Freezing: model.F and model.S already exist, the automated model reformulation will be not performed. \n')
            time_freezing_model_exists = 1;
        end
    else
        if exist('F')
            fprintf('Info on Time-Freezing: model.F provided, model.S missing, the automated model reformulation will be not performed. \n')
            time_freezing_model_exists = 1;
        else
            time_freezing_model_exists = 0;
        end
    end

    %% Check is there a switching function
    if ~exist('f_c')
        if ~exist('c')
            error('Please provide a scalar constraint function model.f_c.')
        else
            f_c = c;
            model.f_c = f_c;
        end
    end
    n_unilateral = length(f_c);

    if n_unilateral > 2
        error('Time-Freezing is currently implemented for only up to two unilateral constraints.')
    end
    %     if n_unilateral > 1
    %         error('Time-Freezing is currently supported only for a single scalar constraint.')
    %     end

    % check dimensions of contacts
    if ~exist('n_dim_contact')
        if  exist('tangent1') && exist('tangent2')
            n_dim_contact = 3;
            fprintf('Info on Time-Freezing: dimension of contact not specified, setting default value n_dim_contact = 3 (since two tangents are provided). \n')
            n_c = 1;
        else
            n_dim_contact = 2;
            fprintf('Info on Time-Freezing: dimension of contact not specified, setting default value n_dim_contact = 2. \n')
        end
    end

    %% Model Checks
    % Check does coeffiecnt of restituion and friction exist and do they have appropiate values, check if f and M exist.
    if ~time_freezing_model_exists
        if exist('e')
            coefficient_of_restitution = e;
        end

        if ~exist('coefficient_of_restitution')
            error('Time-Freezing: Please provide a coefficient of restitution via model.e or model.coefficient_of_restitution. \n')
        end

        if coefficient_of_restitution < 0
            error('Time-Freezing: Please provide a nonnegative model.coefficient_of_restitution or model.e.')
        end

        if abs(1-e)>1 || e<0
            error('The coefficient of restitution e should be in [0,1].')
        end

        if ~exist('mu')
            mu = 0;
        end

        if mu<0
            error('The coefficient of friction mu should be nonnegative.')
        end

        if mu > 0
            friction_is_present = 1;
        else
            friction_is_present = 0;
        end
        % parameter for auxiliary dynamics
        if ~exist('model.a_n') && ~exist('a_n')
            a_n  = 100;
        end
        % no friction for one dimensional models
        if n_dim_contact < 2
            mu = 0;
            friction_is_present  = 0;
            % no friction for a point as there is no tangentional direction
        end

        % dimensions and state space split
        casadi_symbolic_mode = model.x(1).type_name();
        n_x = size(x,1);
        n_q = n_x/2;
        if ~exist('q','var') && ~exist('v','var')
            q = x(1:n_q);
            v = x(n_q+1:2*n_q);
        end
        % update model with this data
        model.n_q = n_q; model.q = q; model.v = v;





        % check model function
        if ~exist('f','var')
            error('The function model.f, in dv/dt =  f(q,v) + ... is not provided to the model.');
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
            fprintf('Info on Time-Freezing: Inertia matrix M is not provided, set to default value: M = eye(n_q). \n')
            M = eye(n_q);
            invM = inv(M);
            model.M = M;
        else
            invM = inv(M);
        end

        %         % is the inertia matrix inverted and treated directly or lifted
        if settings.time_freezing_lift_forces
            invM = eye(n_q); % inverse not used when lifting.
        else
            f = invM*f; % update model by mutpliting it with the inverse
            model.f = model.f;
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



        % multiple impacts
        if n_unilateral > 2
            friction_is_present = 0;
            n_dim_contact = 1;
            % remark: friction is currently only if a single constraint is present.
            % remark : multiple contact are only for planar contacts avilable
        end
        %% Virtual forces
        % add control variables that should help the convergence but
        % penalized to be zero in the solution, or by tihght bounds
        if settings.virtual_forces
            n_u_virtual = n_x;
            u_virtual = define_casadi_symbolic(settings.casadi_symbolic_mode,'u_virtual',n_q); % homotopy parameter;

            if settings.virtual_forces_convex_combination
                psi_vf = define_casadi_symbolic(settings.casadi_symbolic_mode,'psi_vf');
                model.psi_vf = psi_vf ;
            end

            u = [u;u_virtual];
            f_u_virtual = [zeros(n_q,1);u_virtual;zeros(n_quad+1,1)];
            f_q_virtual = u_virtual'*u_virtual;
            f_q_virtual_fun = Function('f_q_virtual_fun',{u},{f_q_virtual});
            f_u_virtual_fun = Function('f_u_virtual_fun',{u},{f_u_virtual});

            if ~settings.virtual_forces_in_every_mode
                if settings.virtual_forces_convex_combination
                    f = (1-psi_vf)*f+psi_vf*(u_virtual+f_gravity);
                else
                    f = f+[u_virtual;zeros(n_quad,1)];
                end
            end
            model.u = u;
            model.f_u_virtual_fun  = f_u_virtual_fun;
            model.f_q_virtual_fun  = f_q_virtual_fun;
        end

        %% Clock state and dimensions
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

        % compute normal vector
        if ~exist('nabla_q_f_c','var')
            nabla_q_f_c = f_c.jacobian(q)';
        end

        if is_zero(nabla_q_f_c)
            error('The normal vector should have at least one non-zero entry.')
        end
        F = {}; c = {};  S = {};

        %% Time-Freezing for mechanical impacts
        if coefficient_of_restitution == 0
            settings.time_freezing_inelastic = 1; % flag tha inealstic time-freezing is using (for hand crafted lifting)
            settings.pss_mode = 'Step'; % currently only with step mode supported
            % unconstrained ODE
            f_ode = [v;f;1];

            % auxiliary dynamics for normal velocity jumps;
            % create auxiliary dynamics
            switch n_unilateral
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
            if friction_is_present
                switch n_dim_contact
                    case 2
                        switch n_unilateral
                            case 1
                                if exist('tangent1','var') || exist('tangent','var')
                                    fprintf('Time-freezing: tangents provided by the user. \n');
                                    if exist('tangent','var')
                                        tangent1 = tangent;
                                    end
                                else
                                    fprintf('Time-freezing: user did not provide tangent vectors at contact, generating them automatically...\n');
                                    t1 = 1; t2 = 1;
                                    if is_zero(nabla_q_f_c(1))
                                        t2 = 0;
                                    elseif is_zero(nabla_q_f_c(2))
                                        t1 = 0;
                                    else
                                        t2 = -nabla_q_f_c(1)*t1/nabla_q_f_c(2);
                                    end
                                    tangent1 = [t1;t2];
                                end
                            case 2
                                if ~exist('tangent1','var') || ~exist('tangent2','var')
                                    error('Time-freezing: user did not provide tangent vectors at contact points. Please specify model.tangent1 and model.tangent2.');
                                end

                        end
                    case 3
                        if exist('tangent1') && exist('tangent2')
                            fprintf('Time-freezing: tangents provided by the user. \n');
                        else
                            % tangent 1
                            %                             fprintf('Time-freezing: user did not provide tangent vectors at contact, generating them automatically...\n');
                            fprintf('Time-freezing: user did not provide tangent vectors at contact, settings mu = 0. \n');
                            mu = 0;
                            t1 = 1; t2 = 1; t3 = 1;
                            ind_tangent1 = 1;
                            if ~is_zero(nabla_q_f_c(1))
                                t1 = -(t2*nabla_q_f_c(2)+t3*nabla_q_f_c(3))/nabla_q_f_c(1);
                                ind_tangent1 = 1;
                            elseif ~is_zero(nabla_q_f_c(2))
                                t2 = -(t1*nabla_q_f_c(1)+t3*nabla_q_f_c(3))/nabla_q_f_c(2);
                                ind_tangent1 = 2;
                            else
                                t3 = -(t1*nabla_q_f_c(1)+t2*nabla_q_f_c(2))/nabla_q_f_c(3);
                                ind_tangent1 = 3;
                            end
                        end
                end

                % create auxiliary dynamics
                a_t = mu*a_n;
                switch n_dim_contact
                    case 2
%                         switch n_u_virtual
%                             case 1
                                f_aux_t1 = [zeros(n_q,1);invM*tangent1*a_t;zeros(n_quad+1,1)];
                                f_aux_t2 = -[zeros(n_q,1);invM*tangent1*a_t;zeros(n_quad+1,1)];
%                             case 2
%                                 f_aux_t1 = [zeros(n_q,1);invM*tangent1*a_t;zeros(n_quad+1,1)];
%                                 f_aux_t2 = -[zeros(n_q,1);invM*tangent1*a_t;zeros(n_quad+1,1)];
%                                 f_aux_t3 = [zeros(n_q,1);invM*tangent2*a_t;zeros(n_quad+1,1)];
%                                 f_aux_t4 = -[zeros(n_q,1);invM*tangent2*a_t;zeros(n_quad+1,1)];
%                         end
                        
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
            if friction_is_present
                switch n_dim_contact
                    case 2
                        c3 = tangent1'*v;
                    case 3
                        c3 = eps_tangential-v_tan_norm;
                end
            else
                c3 = [];
            end

            if ~friction_is_present
                switch n_unilateral
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
                switch n_unilateral
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
            if settings.virtual_forces && settings.virtual_forces_in_every_mode
                if settings.virtual_forces_convex_combination
                    F_temp = F{1};
                    F_temp(n_q+1:2*n_q,:) = (1-psi_vf)*F_temp(n_q+1:2*n_q,:)+psi_vf*(u_virtual+f_gravity);
                    F{1} = F_temp;
                else
                    F{1} = F{1}+f_u_virtual;
                end
            end
            % store results and flag that model is created
            model.F = F; model.c = c;  model.S = S;
            time_freezing_model_exists = 1;
        else
            % elastic
            if ~exist('k_aux')
                k_aux = 10;
                if settings.print_level > 1
                    fprintf('Info on Time-Freezing: Setting default value for k_aux = 10.\n')
                end
            end
            temp1 =2*abs(log(coefficient_of_restitution));
            temp2 = k_aux/(pi^2+log(coefficient_of_restitution)^2);
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
else
    fprintf('Info on Time-Freezing: No action was done. Consider setting settings.time_freezing = 1, if calling this function.\n')
end
settings.time_freezing_model_exists = time_freezing_model_exists;
settings.friction_is_present  = friction_is_present;
model.n_dim_contact = n_dim_contact;
model.n_unilateral = n_unilateral;
end