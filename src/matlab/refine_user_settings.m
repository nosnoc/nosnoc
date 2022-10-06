function [settings] = refine_user_settings(settings)
unfold_struct(settings,'caller')

%% Sanity check of RK Schmes
if ~any(strcmp(list_of_all_rk_schemes,irk_scheme))
    if print_level >=1
        fprintf(['Info: The user provided RK scheme: ' irk_scheme ' is not supported, switching to Radau-IIA. \n'])
        fprintf(['See settings.list_of_all_rk_schemes in the default settings for an overview. \n'])
    end
    irk_scheme = 'Radau-IIA';
    settings.irk_scheme  = irk_scheme ;
end

if any(strcmp(list_of_all_rk_schemes(5:end),irk_scheme))
    if print_level >=1
        fprintf(['Info: The user provided RK scheme: ' irk_scheme ' is only avilabile in the differential representation.\n']);
    end
    irk_representation = 'differential';
    settings.irk_representation = irk_representation;
end

%% corrections for pss_mode
if isequal(settings.pss_mode,'stewart')  || isequal(settings.pss_mode,'stwrt') || isequal(settings.pss_mode,'indicator')
    settings.pss_mode = 'Stewart';
end

if isequal(settings.pss_mode,'step')  || isequal(settings.pss_mode,'stp') || isequal(settings.pss_mode,'Heaviside') ...
        || isequal(settings.pss_mode,'heaviside') || isequal(settings.pss_mode,'AP')
    settings.pss_mode = 'Step';
end

if ~isequal(settings.pss_mode,'Stewart') && ~isequal(settings.pss_mode,'Step')
    error('Please use for settings.pss_mode either ''Stewart''  or  ''Step''.' );
end

if settings.n_depth_step_lifting   < 2
    settings.n_depth_step_lifting  = 2;
        if print_level >=1
            fprintf(['Info: n_depth_step_lifting should be at least 2, changing n_depth_step_lifting  = 2.\n']);
        end
end
%% Some settings refinements
% update print_level
if print_level <4
    settings.opts_ipopt.ipopt.print_level=0;
    settings.opts_ipopt.print_time=0;
    settings.opts_ipopt.ipopt.sb= 'yes';
elseif print_level == 4
    settings.opts_ipopt.ipopt.print_level=0;
    settings.opts_ipopt.print_time=1;
    settings.opts_ipopt.ipopt.sb= 'no';
else
    settings.opts_ipopt.ipopt.print_level = 5;
end

if settings.time_freezing
    settings.local_speed_of_time_variable = 1;
end

%% preprocess settings refined
if settings.h_fixed_change_sigma == 0
   settings.h_fixed_max_iter = 1; 
end
end

