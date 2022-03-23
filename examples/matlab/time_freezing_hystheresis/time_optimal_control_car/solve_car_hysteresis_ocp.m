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
function [varargout] = solve_car_hysteresis_ocp(varargin)
% crate and solve and OCP for an example time freezing problem
%%
import casadi.*
if nargin == 0
    [settings] = default_settings_fesd();  %% Optionally call this function to have an overview of all options.
elseif nargin == 1
    settings = varargin{1};
    unfold_struct(settings,'caller');
    model = [];
else
    settings = varargin{1};
    model = varargin{2};
    unfold_struct(model,'caller');
    unfold_struct(settings,'caller');
end
%% Complete Settings
settings = fill_in_missing_settings(settings);
%% Generate Model
model = car_hystheresis_model_voronoi(model);
%% Formulate NLP;
[solver,solver_initalization, model,settings] = create_nlp_fesd(model,settings);
%% Get variables into main workspace
unfold_struct(model,'caller');
unfold_struct(settings,'caller');
unfold_struct(solver_initalization,'caller');
%% Solve NLP
try
    [sol,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);
    solver_message = sprintf('Test sucessful with complementarity residual = %2.2e. \n',stats.complementarity_stats(end));
    results.status = 1;
        if stats.complementarity_stats(end) > 1e-3
            solver_message = sprintf('Homotopy loop terminated with large complementarity residual: %2.2e . \n',stats.complementarity_stats(end));
            results.status = 0;
        end
catch
    solver_message = 'Test faild. Homotopy loop diverged,\n';
    results.status = 0;
    sol = [];
    stats = [];
end
%% Process results
if results.status == 1
    w_opt = full(sol.x);
    diff_states = w_opt(ind_x);
    for i = 1:n_x
        eval( ['x' num2str(i) '_opt = diff_states(' num2str(i) ':n_x+n_x*n_s:end);']);
    end
    % Check the objective
    if ~isempty(ind_t_final)
        T_opt = w_opt(ind_t_final);
        clocl_state_diff = x5_opt(end) - T_opt;
    else
        T_opt = nan;
    end
    T_phy = x5_opt(end);
    time_message = sprintf('Phyisical time T_phy =  %2.3f s, numerical time T_num = %2.3f s, final time (if any) T_final = %2.3f s.\n',T_phy,T,T_opt);
    terminal_error = sprintf('Terminal constraint satisfication:= %2.3e .\n',norm([q_goal;v_goal]-[x1_opt(end);x2_opt(end)]));
    messages.time_message = time_message;
    messages.terminal_error = terminal_error;
end
%% Output
messages.solver_message = solver_message;
if results.status
    results.sol = sol;
    results.w_opt = full(sol.x);
    if store_all_homotopy_iterates
        results.W = sol.W;
    end
end
results.messages = messages;
% varargout
varargout{1} = results;
varargout{2} = stats;
varargout{3} = solver_initalization;
varargout{4} = settings;
varargout{5} = model;
end

