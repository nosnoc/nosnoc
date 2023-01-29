% BSD 2-Clause License

% Copyright (c) 2023, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

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
classdef NosnocFormulationObject < handle
% NosnocFormulationObject
    properties
        % Primal variables
        w
        w0
        lbw
        ubw
        % Constraints
        g
        lbg
        ubg
        % Cost
        cost
        % Objective
        objective
    end

    methods
        function obj = NosnocFormulationObject()
            import casadi.*

            % Primal variables.
            obj.w = [];
            obj.w0 = [];
            obj.lbw = [];
            obj.ubw = [];

            % Constraints
            obj.g = [];
            obj.lbg = [];
            obj.ubg = [];

            % Cost
            obj.cost = 0;

            % Objective
            obj.objective = 0;
        end

        function [] = addVariable(obj, symbolic, type, lb, ub, initial, varargin)

            lens = [size(symbolic,1),size(lb,1),size(ub,1), size(initial,1)];
            if ~all(lens == lens(1))
                symbolic
                lb
                ub
                initial
                error("mismatched dims")
            end
            
            n = size(symbolic, 1);
            n_w = length(obj.w);

            obj.w = vertcat(obj.w, symbolic);
            obj.lbw = [obj.lbw; lb];
            obj.ubw = [obj.ubw; ub];
            obj.w0 = [obj.w0; initial];

            new_indices = (n_w+1):(n_w+n);

            if length(varargin) ~= 0
                obj.(strcat('ind_', type)){varargin{:}} = new_indices;
            else
                obj.(strcat('ind_', type)) = [obj.(strcat('ind_', type)),new_indices];
            end
        end

        function [] = addConstraint(obj, symbolic, varargin)
            import casadi.*

            p = inputParser;
            % TODO: add checks.
            addRequired(p, 'obj');
            addRequired(p, 'symbolic');
            addOptional(p, 'lb', []);
            addOptional(p, 'ub', []);
            addOptional(p, 'type', []);
            parse(p, obj, symbolic, varargin{:});
            
            n = length(symbolic);
            n_g = length(obj.g);

            if ismember('lb', p.UsingDefaults)
                lb = zeros(n, 1);
            else
                lb = p.Results.lb;
            end
            if ismember('ub', p.UsingDefaults)
                ub = zeros(n, 1);
            else
                ub = p.Results.ub;
            end

            obj.g = vertcat(obj.g, symbolic);
            obj.lbg = [obj.lbg; lb];
            obj.ubg = [obj.ubg; ub];

            if ~ismember('type', p.UsingDefaults)
                new_indices = (n_g+1):(n_g+n);

                obj.(strcat('ind_', type)) = [obj.(strcat('ind_', type)), new_indices];
            end
        end
    end
end
