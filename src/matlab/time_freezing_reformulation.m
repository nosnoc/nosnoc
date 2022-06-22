function [model,settings] = time_freezing_reformulation(model,settings)
import casadi.*
%% Load settings and model details
unfold_struct(model,'caller');
time_freezing = settings.time_freezing;
time_freezing_model_exists = 0;
if time_freezing
    % sanity heck
    %% check is the model partially filled
    if exist('S')
        if exist('F')
            fprintf('Info on Time-Freezing: model.F and model.S already exist, the automated model reformulation will be not started. \n')
            time_freezing_model_exists = 1;
        end
    else
        if exist('F')
            fprintf('Info on Time-Freezing: model.F provided, model.S missing, the automated model reformulation will be not started. \n')
            time_freezing_model_exists = 1;
        else
            time_freezing_model_exists = 0;
        end
    end
    % CoR check
    %% Check is there a switching function
    if ~exist('f_c')
        if ~exist('c')
            error('Please provide a scalar constraint function model.f_c.')
        else
            f_c = c;
            model.f_c = f_c;
        end
    end
    if length(f_c) > 1
        error('Time-Freezing is currently supported only for a single scalar constraint.')
    end
    % check dimensions of contacts
    if ~exist('n_dim_contact')

        if  exist('tangent1') && exist('tangent2')
            n_dim_contact = 3;
            fprintf('Info on Time-Freezing: dimension of contact not specified, setting default value n_dim_contact = 3 (since two tangents are provided). \n')
        else
            n_dim_contact = 2;
            fprintf('Info on Time-Freezing: dimension of contact not specified, setting default value n_dim_contact = 2. \n')
        end
    end

    %% Check does coeffiecnt of restituion and friction exist and do they have appropiate values
    if ~time_freezing_model_exists
        if exist('e')
            coefficient_of_restitution = e;
        end

        if ~exist('coefficient_of_restitution')
            error('Time-Freezing: Please provide a coefficient of restitution via model.e or model.coefficient_of_restitution. \n')
        end

        if coefficient_of_restitution < 0
            error('Time-Freezing: Please provide a positive model.coefficient_of_restitution or model.e.')
        end

        if abs(1-e)>1 || e<0
            error('The coefficient of restitution e should be in [0,1].')
        end
        if ~exist('mu')
            mu = 0;
        end
        if  mu<0
            error('The coefficient of friction mu should be non-negative.')
        end

        if mu > 0
            friction_is_present = 1;
        else
            friction_is_present = 0;
        end
        %% Dimensions, states and clock state
        casadi_symbolic_mode = model.x(1).type_name();
        t = define_casadi_symbolic(casadi_symbolic_mode,'t',1);

        % update lower and upper bounds of lbx and ubx
        if exist('lbx')
            model.lbx = [model.lbx;-inf];
        end
        if exist('ubx')
            model.ubx = [model.ubx;inf];
        end

        n_x = size(x,1);
        n_q = n_x/2;
        if ~exist('q','var') && ~exist('v','var')
            q = x(1:n_q);
            v = x(n_q +1:end);
        end



        x = [x;t];
        model.x = x;
        model.x0 = [model.x0;0];

        % check model function
        if ~exist('f','var')
            error('The function model.f, in dv/dt =  f(q,v) + ... is not provided to the model.');
        end

        % Check intertia matrix

        if ~exist('M','var')
            fprintf('Info on Time-Freezing: Neither M or G are provided, set to default value: M = eye(n_q). \n')
            M = eye(n_q);
            model.M = M;
        else
            invM = inv(M);
        end
        % compute normal vector
        if ~exist('nabla_q_f_c','var')
            nabla_q_f_c = f_c.jacobian(q)';
        end
        if is_zero(nabla_q_f_c)
            error('The normal vector should have at least one non-zero entry.')
        end

        F = {};
        c = {};
        S = {};

        %% Time-Freezing for mechanical impacts
        if coefficient_of_restitution > 0
            % elastic
            if ~exist('k_aux');
                k_aux = 10;
                if  settings.print_level > 1;
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
            f_aux_n = N*K*N'*[q;v];

            f = [f;1];
            f_aux_n = [f_aux_n;0];

            % updated with clock state
            model.f = f;
            model.F = [f, f_aux_n];
            model.S = [1; -1];

            time_freezing_model_exists = 1;
        else
            % inelastic
            settings.time_freezing_inelastic = 1; % flag tha inealstic time-freezing is using (for hand crafted lifting)
            settings.pss_mode = 'Step'; % currently only with step mode supported
            % compute tangents
            if n_dim_contact< 2
                mu = 0;
                friction_is_present  = 0;
                % no friction for a point as there is no tangentional direction
            end
            if friction_is_present  && n_dim_contact > 1
                switch n_dim_contact
                    case 2
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
                    case 3
                        if exist('tangent1') && exist('tangent2')
                            fprintf('Time-freezing: tangents provided by the user. \n');
                        else
                            % tangent 1
                            fprintf('Time-freezing: user did not provide tangent vectors at contact, generating them automatically...\n');
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
                            %     tangent1 = [t1;t2;t3];
                            %     b = ones(3,1);
                            %     b(ind_tangent1) = -
                        end
                end
            end

            %% create auxiliary dynamics
            if ~exist('model.a_n') && ~exist('a_n')
                a_n  = 100;
            end
            f_aux_n = [nabla_q_f_c zeros(n_q,1);...
                zeros(n_q,1) invM*nabla_q_f_c]*[0;a_n];
            if friction_is_present
                a_t = mu*a_n;
                f_aux_t1 = [tangent1 zeros(n_q,1);...
                    zeros(n_q,1) invM*tangent1]*[0;a_t];
                if n_dim_contact> 2
                    f_aux_t2 = [tangent2 zeros(n_q,1);...
                        zeros(n_q,1) invM*tangent2]*[0;a_t];
                end
            end

            c1 = f_c;
            c2 = nabla_q_f_c'*v;

            if friction_is_present
                switch n_dim_contact
                    case 2
                        c3 = tangent1'*v;
                    case 3
                        c3 = tangent1'*v;
                        c4 = tangent2'*v;
                end
            end

            % PSS ode formulation
            f_ode = [v;f;1];
            f_aux_n = [f_aux_n;0];

            F{1} = [f_ode f_aux_n];
            c{1} = [c1;c2];
            S{1} = [1 0;-1 1;-1 -1];
            % frictional impact
            if friction_is_present
                switch n_dim_contact
                    case 2
                        f_aux_t1 = [f_aux_t1;0];
                        F{1} = [f_ode ...
                                f_aux_n+f_aux_t1 ...
                                f_aux_n-f_aux_t1];
                        c{1} = [c1;c2;c3];
                        S{1} = [1 0 0;-1 -1 -1;-1 -1 1];
                    case 3
                        f_aux_t1 = [f_aux_t1;0];
                        f_aux_t2 = [f_aux_t2;0];
                        F{1} = [f_ode ...
                                f_aux_n+f_aux_t1+f_aux_t2...
                                f_aux_n+f_aux_t1-f_aux_t2...
                                f_aux_n-f_aux_t1+f_aux_t2...
                                f_aux_n-f_aux_t1-f_aux_t2 ];
                        c{1} = [c1;c2;c3;c4];
                        S{1} = [1 0 0 0;...
                            -1 -1 -1 -1;...
                            -1 -1 -1 1;...
                            -1 -1 1 -1;...
                            -1 -1 1 1];
                end
            end
            model.F = F;
            model.c = c;
            model.S = S;
            time_freezing_model_exists = 1;
        end
    end
else
    fprintf('Info on Time-Freezing: No action was done. Consider setting settings.time_freezing = 1, if calling this function.\n')
end
settings.time_freezing_model_exists = time_freezing_model_exists;
settings.friction_is_present  = friction_is_present;
model.n_dim_contact = n_dim_contact;
end