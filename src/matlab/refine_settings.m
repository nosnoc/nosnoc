function [settings] = refine_settings(settings);

%% Unfold user structure
unfold_struct(settings,'caller')

%% MPCC and Homotopy
% sigma_k = sigma_0;

% Modify MPCC Settings
if mpcc_mode == 1
   N_homotopy = 1;
   sigma_0 = 0;
   sigma_N = 0;
   mpcc_mode = 3;
end

if mpcc_mode == 4 || mpcc_mode == 10
%    pointwise_or_integral = 0;
   cross_complementarity_mode = 10;
end

if mpcc_mode >= 8 && mpcc_mode <= 10;
%     s_elastic_min = 1e-16;
    s_elastic_max = inf;
    N_homotopy = 1;   
end

if equidistant_control_grid == 0
    couple_across_stages = 1;
end

if nonlinear_sigma_rho_constraint
    convex_sigma_rho_constraint = 1;
end

%% Correct contradictring settings, complete missing data (if any)
step_equilibration = step_equilibration*use_fesd;

there_exist_free_x0 = exist('ind_free_x0');
%
% if N_finite_elements == 1
%     equidistant_control_grid = 0;
%     warning('N_finite_elements = 1 and equidistant_control_grid = 1 lead to a equidistant discretization grid. Problem might be infeasible.\n Setting equidistant_control_grid = 1.')
% end
% This option below makes only sense if the h_{n,m} are variable and this is the se if FESD is used.

%% Time Scaling

if (time_freezing || time_optimal_problem) == 1
    time_rescaling = 1;
else
    if time_rescaling
        warning('Time rescaling makes only sense if either time freezing or time optimal problem is on. Setting time_rescaling =0.')
        time_rescaling = 0;
    else
        time_rescaling = 0;
    end
end

if time_rescaling == 0
    use_speed_of_time_variables  = 0;
end

if use_speed_of_time_variables == 0
    local_speed_of_time_variable = 0;
end

%% Impose time
% this gets active only if time freezing is activated and the user had provided impose_terminal_phyisical_time = 0.
if (time_freezing && ~impose_terminal_phyisical_time) 
    time_rescaling = 0;
    % But what if I want to get time_optimal_problem and minimize the
    % numericla time? The implementation of this should be avoided.
end
if impose_terminal_phyisical_time == 0
    warning ('impose_terminal_phyisical_time = 0 is not recommended. It means T \neq T_phy (or T_final \neq T_phy). It is only supported for nonequdistant control grids \n')
end

%% N_finite_elements per stage
%% Number of finite elements
% make a vector
% if length(settings.N_finite_elements) == 1  
%     settings.N_finite_elements = settings.N_finite_elements*ones(settings.N_stages,1)
% elseif length(settings.N_finite_elements) > 1 && length(settings.N_finite_elements) < settings.N_stages
%     settings.N_finite_elements = settings.N_finite_elements(:); % make sure it is a column vector
%     settings.N_finite_elements = [settings.N_finite_elements;settings.N_finite_elements(end)*ones(settings.N_stages-length(settings.N_finite_elements),1)];
% end
if length(N_finite_elements) == 1  
    N_finite_elements = N_finite_elements*ones(N_stages,1);
elseif length(N_finite_elements) > 1 && length(N_finite_elements) < N_stages
    N_finite_elements = N_finite_elements(:); % make sure it is a column vector
    N_finite_elements = [N_finite_elements;N_finite_elements(end)*ones(N_stages-length(N_finite_elements),1)];
end




%% Save data for output into struct
% settings = [];
names = who;
for ii = 1:length(names)
    eval([ 'settings.' names{ii} '=' names{ii} ';'])
end
end