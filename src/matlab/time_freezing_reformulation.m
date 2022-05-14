function [model,settings] = time_freezing_reformulation(model,settings)
import casadi.*
%% Load settings and model details

% model = check_is_user_model_valid(model,settings);
unfold_struct(model,'caller');
time_freezing = settings.time_freezing;
time_freezing_model_exists = 0;
if time_freezing
    % sanity heck
    %% check is the model partially filled

    if exist('S')
        if exist('F')
            fprintf('Info on Time-Freezing: model.F and model.S already exist, the automated model reformulation will be not started. \n.')
            time_freezing_model_exists = 1;
        end
    else
        if exist('F')
            fprintf('Info on Time-Freezing: model.F provided, model.S missing, the automated model reformulation will be not started. \n.')
            time_freezing_model_exists = 1;
        else
            time_freezing_model_exists = 0;
        end
    end
    % CoR check


    %% test are paseed, starting reformulation
    if ~time_freezing_model_exists
        if exist('e')
            coefficient_of_restitution = e;
        end

        if ~exist('coefficient_of_restitution')
            error('Time-Freezing: Please provide a coefficient_of_restitution.')
        end

        if coefficient_of_restitution < 0
            error('Time-Freezing: Please provide a positive coefficient_of_restitution.')
        end

        if coefficient_of_restitution > 0
            % elastic
            if ~exist('k_aux');
                k_aux = 10;
                if  settings.print_level > 1;
                    fprintf('Info on Time-Freezing: Setting default value for k_aux = 10.\n')
                end
            end
            % clock state
            eval(['t = ' settings.casadi_symbolic_mode '.sym(''t'');'])
            n_x = size(x,1);
            q = x(1:n_x/2);
            v = x(n_x/2+1:end);

            % todo: apply the more general formuala with derivatives
            %             nabla_c = c.jacobian(x(1:end-1))*x;
%             c_aux = 2*abs(log(coefficient_of_restitution))/sqrt(k_aux/(log(coefficient_of_restitution)^2+pi^2));
            temp1 =2*abs(log(coefficient_of_restitution));
            temp2 = k_aux/(pi^2+log(coefficient_of_restitution)^2);
            c_aux = temp1/sqrt(temp2);
%                 c_aux = 0.211989;
            f_aux = [0;v(2);0;-k_aux*q(2)-c_aux*v(2);0];
            x = [x;t];
            f = [f;1];
            model.F = [f, f_aux];
            model.S = [1; -1];
            % updated with clock state
            model.f = f;
            model.x = x;
            model.x0 = [model.x0;0];
            time_freezing_model_exists = 1;
        else
            % inelastic
            error('work in progres ...');
        end
    end
else
    fprintf('Info on Time-Freezing: No action was done. Consider setting settings.time_freezing = 1, if calling this function.\n')
end
settings.time_freezing_model_exists = time_freezing_model_exists;

end

