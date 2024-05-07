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

% TODO prehaps it is necessary to have a separate NLP relaxation and direct mpcc solver class.
classdef NosnocSolver < handle
    properties
        mpcc     % Nosnoc Mpcc
        solver_options  % Nosnoc Options
        solver % Nosnoc solver function
    end

    methods
        function obj = NosnocSolver(mpcc, solver_options, solver_type)
            import casadi.*
            if ~exist("solver_type")
                solver_type = 'scholtes_ineq';
            end
            obj.mpcc = mpcc;
            obj.solver_options = solver_options;
            solver_options.assume_lower_bounds = true;
            opts = solver_options;
            mpcc_struct = struct;
            mpcc_struct.x = mpcc.w;
            mpcc_struct.g = mpcc.g;
            mpcc_struct.p = mpcc.p;
            mpcc_struct.f = mpcc.augmented_objective;
            mpcc_struct.G = mpcc.cross_comps(:,1);
            mpcc_struct.H = mpcc.cross_comps(:,2);
            
            obj.solver = nosnoc.solver.mpccsol('nosnoc_solver', solver_type, mpcc_struct, opts);
        end

        function set(obj, type, val)
            mpcc = obj.mpcc;
            if strcmp(type, 'w0')
                if length(obj.mpcc.w0) == length(val)
                    obj.mpcc.w0 = val;
                else
                    error("nosnoc: if initializing w0 all at once you need to provide a vector of corresponding size.")
                end
            else
                % This line uses the index sets collected during the creation of the NosnocProblem and automatically gets
                % the one of the form 'ind_<type>'. This makes this set generic for all variable types.
                ind = obj.mpcc.(strcat('ind_', type));
                if iscell(ind)
                    flat_ind = sort([ind{:}]);
                else
                    flat_ind = ind;
                end
                % TODO: Interpolation

                if iscell(val)
                    % If the passed value is an N_stage by 1 cell array we assume this initialization is done stage wise
                    if ismatrix(val) && size(val, 1) == mpcc.problem_options.N_stages && size(val,2) == 1
                        for ii=1:mpcc.problem_options.N_stages
                            % All variables of each stage are set to the same value
                            for v=ind(ii,:,:)
                                % NOTE: isempty check is needed for possibly unused rk-stage level cells (like in the case of rbp, etc.)
                                if ~isempty(v) && length(v{1}) == length(val{ii})
                                    obj.mpcc.w0(v{1}) = val{ii};
                                end
                            end
                        end
                    % Otherwise if we have an initialization of the form N_stages-by-N_fe we do the same but finite-element-wise
                    elseif ismatrix(val) && size(val, 1) == mpcc.problem_options.N_stages && size(val, 2) == mpcc.problem_options.N_finite_elements
                        for ii=1:mpcc.problem_options.N_stages
                            for jj=1:mpcc.problem_options.N_finite_elements
                                % All variables of each finite element are set to the same value
                                for v=ind(ii,jj,:)
                                    % NOTE: isempty check is needed for possibly unused rk-stage level cells (like in the case of rbp, cls, etc.)
                                    if ~isempty(v) && length(v{1}) == length(val{ii,jj})
                                        obj.mpcc.w0(v{1}) = val{ii,jj};
                                    end
                                end
                            end
                        end
                    else
                        error("nosnoc: Initialization either needs to be Stagewise or Finite Element wise, i.e. provide either an N_stage-by-1 or N_stage-by-N_fe cell array");
                    end
                    % Otherwise we assume that we are initializing via a flat array and we simply check for the same length
                else
                    if ismatrix(val) && length(val) == length(flat_ind)
                        obj.mpcc.w0(flat_ind) = val;
                    else
                        error('nosnoc: set should be a cell array or a flat array')
                    end
                end
            end
        end

        function [results,stats] = solve(obj)
            mpcc = obj.mpcc;
            p0 = mpcc.p0; % TODO(@anton) correct?
            results = obj.solver('x0', mpcc.w0,...
                'lbx', mpcc.lbw,...
                'ubx', mpcc.ubw,...
                'lbg', mpcc.lbg,...
                'ubg', mpcc.ubg,...
                'p', p0);
            stats = obj.solver.stats;
            results = obj.extract_results(results);
        end
       
    end

    methods(Access=private)

        % TODO(@anton) most of this will die and be replaced with vdx :)
        function results = extract_results(obj, results)
            import casadi.*
            mpcc = obj.mpcc;
            % Store differential states
            w_opt = full(results.x);
            results.w = w_opt;
            w_mpcc = w_opt;

            names = get_result_names_from_settings(obj.mpcc.problem_options);
            % populate outputs
            for name=names
                results = form_structured_output(mpcc, w_mpcc, name, results);
            end

            % handle x0 properly
            x0 = w_mpcc(mpcc.ind_x0);
            results.x = [x0, results.x];
            results.extended.x = [x0, results.extended.x];

            if ~isempty(mpcc.ind_v_global)
                v_global = w_mpcc(mpcc.ind_v_global);
                results.v_global = v_global;
            end
            
            u = w_mpcc([mpcc.ind_u{:}]);
            u = reshape(u,obj.mpcc.model.dims.n_u,obj.mpcc.problem_options.N_stages);

            results.u = u;

            if mpcc.problem_options.time_optimal_problem
                T_opt = w_mpcc(mpcc.ind_t_final);
            else
                T_opt = mpcc.problem_options.T;
            end
            results.T = T_opt;

            if obj.mpcc.problem_options.use_fesd
                h_opt = w_mpcc(flatten_ind(mpcc.ind_h))';
            else
                h_opt = [];
                % if settings.time_optimal_problem && ~settings.use_speed_of_time_variables
                %     T = T_opt;
                % end
                for ii = 1:obj.mpcc.problem_options.N_stages
                    h_opt = [h_opt,obj.mpcc.problem_options.T/(obj.mpcc.problem_options.N_stages*obj.mpcc.problem_options.N_finite_elements(ii))*ones(1, obj.mpcc.problem_options.N_finite_elements(ii))];
                end
            end
            results.h = h_opt;

            t_grid = cumsum([0,h_opt]);

            if obj.mpcc.problem_options.use_speed_of_time_variables
                s_sot = w_mpcc(flatten_ind(mpcc.ind_sot));
                if ~obj.mpcc.problem_options.local_speed_of_time_variable
                    s_sot = s_sot*ones(obj.mpcc.problem_options.N_stages,1);
                end
                results.s_sot = s_sot;
            end
            
            %% Adapt the grid in case of time optimal problems
            if obj.mpcc.problem_options.time_optimal_problem
                if obj.mpcc.problem_options.use_speed_of_time_variables
                    s_sot = w_mpcc(flatten_ind(mpcc.ind_sot));
                    if ~obj.mpcc.problem_options.local_speed_of_time_variable
                        s_sot = s_sot*ones(obj.mpcc.problem_options.N_stages,1);
                    end
                    results.s_sot = s_sot;
                    h_rescaled = [];
                    ind_prev = 1;
                    for ii = 1:obj.mpcc.problem_options.N_stages
                        h_rescaled = [h_rescaled,h_opt(ind_prev:obj.mpcc.problem_options.N_finite_elements(ii)+ind_prev-1).*s_sot(ii)];
                        ind_prev = ind_prev+obj.mpcc.problem_options.N_finite_elements(ii);
                    end
                    t_grid = cumsum([0,h_rescaled]);
                else
                    t_grid = cumsum([0,h_opt]);
                end
            end
            ind_t_grid_u = cumsum([1; obj.mpcc.problem_options.N_finite_elements]);

            if obj.mpcc.problem_options.dcs_mode == DcsMode.CLS
                x_with_impulse = x0;
                t_with_impulse = kron(t_grid, ones(2,1));
                for ii=1:size(results.structured.x,1)
                    for jj=1:size(results.structured.x,2)
                        x_with_impulse = [x_with_impulse,results.structured.x_left_bp{ii,jj}];
                        x_fe = results.structured.x{ii,jj};
                        x_with_impulse = [x_with_impulse, x_fe(:,end)];
                    end
                end
                results.x_with_impulse = x_with_impulse;
                results.t_with_impulse = t_with_impulse(1:end-1);
            end


            results.t_grid = t_grid;
            results.t_grid_u = t_grid(ind_t_grid_u);

            results.f = full(results.f);
            results.g = full(results.g);
            results.objective = results.f;
        end
    end
end
