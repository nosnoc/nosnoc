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

classdef Snopt < handle % TODO maybe handle not necessary, revisit.
    properties

    end

    methods
        function solver = construct_solver(obj, nlp, solver_options, time_remaining)
            opts_casadi_nlp = solver_options.opts_casadi_nlp;
            opts_casadi_nlp = rmfield(opts_casadi_nlp, 'ipopt');
            opts_casadi_nlp = rmfield(opts_casadi_nlp, 'worhp');
            opts_casadi_nlp = rmfield(opts_casadi_nlp, 'uno');
            if solver_options.timeout_wall
                if exist('time_remaining')
                    opts_casadi_nlp.snopt.Time_limit = time_remaining;
                    opts_casadi_nlp.snopt.Timing_level = 3;
                else
                    opts_casadi_nlp.snopt.Time_limit = solver_options.timeout_wall;
                    opts_casadi_nlp.snopt.Timing_level = 3;
                end
            end

            nlp.create_solver(opts_casadi_nlp, 'snopt');
        end

        function solver_stats = cleanup_solver_stats(obj, solver_stats)
            if ~isfield(solver_stats, 'iterations')
                solver_stats.iterations = [];
            end
        end

        function failed = check_iteration_failed(obj, stats)
            switch stats.solver_stats(end).return_status
                case {'Finished successfully'}
                    failed = false;
                otherwise
                    failed = true;
            end
        end

        function timeout = check_timeout(obj, stats)
            switch stats.solver_stats(end).return_status
                case {'Resource limit error'}
                    if strcmp(stats.solver_stats(end).secondary_return_status, 'time limit reached');
                        timeout = 1;
                    else
                        timeout = 0;
                    end
                otherwise
                    timeout = 0;
            end
        end

        function w_opt = w_opt_from_results(obj, nlp_results)
            w_opt = full(nlp_results.x);
        end
        function f = f_from_results(obj, nlp_results)
            f = full(nlp_results.f);
        end
        function g = g_from_results(obj, nlp_results)
            g = full(nlp_results.g);
        end

        function print_nlp_iter_header(obj)
            fprintf('\n|%-5s|%-10s|%-10s|%-10s|%-10s|%-10s|%-20s|%-30s\n', 'iter', 'sigma', 'compl_res', 'objective', 'CPU time', 'NLP iter', 'status', 'secondary status');
        end
        
        function print_nlp_iter_info(obj, stats)
            solver_stats = stats.solver_stats(end);
            ii = size(stats.solver_stats, 2);

            fprintf('\n|%-5d|%6.2e|%6.2e|%6.2e|%6.3f|%-20s|%-30s',...
                    ii, stats.sigma_k(end), stats.complementarity_stats(end), ...
                    stats.objective(end), stats.cpu_time(end), solver_stats.return_status, solver_stats.secondary_return_status);
        end
    end
end
