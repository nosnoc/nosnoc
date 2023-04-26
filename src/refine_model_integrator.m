% BSD 2-Clause License

% Copyright (c) 2022, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% This file is part of NOSNOC.

function [model] = refine_model_integrator(model,settings)
% This function changes possible inconsistent discretization settings.
%% Unfold user structure
unfold_struct(model,'caller')
%% If different names are used...
if exist('N_stg','var')
    N_stages = N_stg;
    model.N_stage = N_stages;
end

if exist('N_FE','var') 
    N_finite_elements = N_FE;
    model.N_finite_elements  = N_finite_elements ;
end

if exist('N_fe','var') 
    N_finite_elements = N_fe;
    model.N_finite_elements = N_fe;
end



%%
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

if exist("N_sim")
    T = T_sim/N_sim;
    h_sim = T_sim/(N_sim*N_stages*N_finite_elements);
    if settings.print_level >= 2 && exist("h_sim")
        fprintf('Info: N_sim is given, so the h_sim provided by the user is overwritten.\n')
    end
else
    error('Provide N_sim for the integration.')
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