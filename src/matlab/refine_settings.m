% Copyright 2022 Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% This file is part of NOSNOC.

% The 2-Clause BSD License

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
function [settings] = refine_settings(settings)

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
   cross_comp_mode = 12;
end

if mpcc_mode >= 8 && mpcc_mode <= 10;
    s_elastic_max = inf;
    N_homotopy = 1;   
end

if s_elastic_0 > s_elastic_max
    s_elastic_max = s_elastic_0;
    if print_level >= 2
        fprintf('Info: s_elastic_0 > s_elastic_max , setting s_elastic_max = s_elastic_0. \n');
    end
end

if nonlinear_sigma_rho_constraint
    convex_sigma_rho_constraint = 1;
end

%% Handling constraint evaluation. 
if ~exist('x_box_at_fe')
    x_box_at_fe = 1;
end
if ~exist('x_box_at_stg')
    x_box_at_stg = 1;
    if isequal(irk_scheme,'differential') &&  ~lift_irk_differential
        x_box_at_stg = 0;
    end
end

%% Correct contradictring settings, complete missing data (if any)
there_exist_free_x0 = exist('ind_free_x0');

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

if exist('cross_complementarity_mode')
    cross_comp_mode = cross_complementarity_mode;
end



%% Impose time
% this gets active only if time freezing is activated and the user had provided impose_terminal_phyisical_time = 0.
% 
% if time_freezing || stagewise_clock_constraint
%     % using them is the only way to rescale the time;
%     use_speed_of_time_variables = 1;
%     local_speed_of_time_variable = 1;
% end
% 
% if stagewise_clock_constraint
%     time_freezing = 1;
% end

if (time_freezing && ~impose_terminal_phyisical_time) 
    time_rescaling = 0;
end

if impose_terminal_phyisical_time == 0
    warning ('impose_terminal_phyisical_time = 0 is not recommended. It means T \neq T_phy (or T_final \neq T_phy). It is only supported for nonequidistant control grids.')
end

% lifting does not make sense in integral mode
if isequal(irk_representation,'integral') 
    lift_irk_differential = 0;
end

%% Save data for output into struct
% settings = [];
names = who;
for ii = 1:length(names)
    eval([ 'settings.' names{ii} '=' names{ii} ';'])
end
end