%
%    This file is part of NOS-NOC.
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
function [model] = refine_model_integrator(model,settings)
% This function changes possible inconsistent discretization settings.
%% Unfold user structure
unfold_struct(model,'caller')
if ~exist("T_sim")
    error('Please provide the total simulation time in model.T_sim')
end

if ~exist("N_stages")
    if settings.print_level >= 2
        fprintf('Info: number of stages N_stages not provided, setting to defaul value N_stages = 1.\n')
    end
    N_stages = 1;
end


if ~exist("N_finite_elements")
    if settings.print_level >= 2
        fprintf('Info: number of stages N_finite_elements not provided, setting to defaul value N_finite_elements = 2.\n')
    end
    N_finite_elements= 2;
end

if length(N_finite_elements)>1
    if settings.print_level >= 2
        fprintf('Info: Different number of finite elements per stage is not supported in integration, will use the first element of N_finite_elements.\n')
    end
    N_finite_elements = N_finite_elements(1);
end


if N_finite_elements == 1 && N_stages == 1 && settings.use_fesd == 1
    if settings.print_level >= 2
        fprintf('Info: N_finite_elements = 1 and N_stages =1 with using FESD is infeasible, resetting N_finite_elements to 2 .\n')
    end
    N_finite_elements = 2;
end



additional_residual_ingeration_step = 0;
if exist("N_sim")
    T = T_sim/N_sim;
    h_sim = T_sim/(N_sim*N_stages*N_finite_elements);
    if settings.print_level >= 2 && exist("h_sim")
        fprintf('Info: N_sim is given, so the h_sim provided by the user is overwritten.\n')
    end
elseif exist("h_sim")
    N_sim_fractional = T_sim/(h_sim*N_finite_elements*N_stages);
    N_sim = floor(N_sim_fractional);
    T = h_sim*N_finite_elements*N_stages;
    if N_sim_fractional  - N_sim > 1e-16
        additional_residual_ingeration_step = 1;
        T_residual = T_sim-N_sim*h_sim*N_finite_elements*N_stages;
        if settings.print_level >= 2
            fprintf('Info: N_sim = T_sim/h_sim is not an integerd, adding a smaller integration step at the end to reach T_sim.\n')
        end
    end
else
    error('Provide either h_sim or N_sim for the integration.')
end

if N_stages > 3 || N_finite_elements > 3
    if settings.print_level >=2
        fprintf('Info: The user provided N_stages = %d and N_finite_elements = %d, smaller values might lead to faster computations.\n',N_stages,N_finite_elements)
    end
end


settings.equidistant_control_grid = 0;
%% Save data for output into struct
% settings = [];
names = who;
for ii = 1:length(names)
    eval([ 'model.' names{ii} '=' names{ii} ';'])
end

end