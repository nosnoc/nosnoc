%
%    This file is part of NOSNOC.
%
%    NOS-NOC -- A software for NOnSmooth Numerical Optimal Control.
%    Copyright (C) 2022 Armin Nurkanovic, Moritz Diehl (ALU Freiburg).
%
%    NOS-NOC is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 3 of the License, or (at your option) any later version.
%
%    NOS-NOC is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with NOS-NOC; if not, write to the Free Software Foundation,
%    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%
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


if equidistant_control_grid == 0
    couple_across_stages = 1;
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
step_equilibration = step_equilibration*use_fesd;
there_exist_free_x0 = exist('ind_free_x0');

    % if use_fesd & N_finite_elements < 2
    %     N_finite_elements = 2;
    %     warrning('Info:In FESD there must be at least 2 finite elements per control interval. Setting N_finite_elements = 2.')
    % end
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
    warning ('impose_terminal_phyisical_time = 0 is not recommended. It means T \neq T_phy (or T_final \neq T_phy). It is only supported for nonequdistant control grids \n')
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