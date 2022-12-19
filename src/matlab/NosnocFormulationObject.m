classdef (Abstract) NosnocFormulationObject < handle
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
    end

    methods
        function obj = NosnocFormulationObject()
            import casadi.*

            % Primal variables.
            obj.w = SX([]);
            obj.w0 = [];
            obj.lbw = [];
            obj.ubw = [];

            % Constraints
            obj.g = SX([]);
            obj.lbg = [];
            obj.ubg = [];

            % Cost
            obj.cost = SX.zeros(1);
        end

        function [] = addVariable(obj, symbolic, idx, lb, ub, initial, varargin)
            import casadi.*

            p = inputParser();
            p.FunctionName('addVariable')
            
            % TODO: add checks.
            addRequired(p, 'obj');
            addRequired(p, 'symbolic');
            addRequired(p, 'lb');
            addRequired(p, 'ub');
            addRequired(p, 'initial');
            addRequired(p, 'idx');
            addOptional(p, 'stage', []);
            addOptional(p, 'sys', []);
            parse(p, obj, symbolic, lb, ub, initial, idx, varargin{:});
            
            n = length(symbolic);
            nw = length(obj.w);

            obj.w = vertcat(obj.w, symbolic);
            obj.lbw = [obj.lbw; lb];
            obj.ubw = [obj.ubw; ub];
            obj.w0 = [obj.w0; initial];

            new_indices = (nw+1):(nw+n);
            
            if ismember('stage',p.UsingDefaults)
                obj.(strcat('ind_', idx)) = [obj.(strcat('ind_', idx)), new_indices];
            else
                if ~ismember('sys',p.UsingDefaults)
                    obj.(strcat('ind_', idx)){stage, sys} = new_indices;
                else
                    obj.(strcat('ind_', idx)){stage} = new_indices;
                end
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
            addOptional(p, 'idx', []);
            parse(p, obj, symbolic, varargin{:});
            
            n = length(symbolic);
            ng = length(obj.g);

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

            if ~ismember('idx', p.UsingDefaults)
                new_indices = (ng+1):(ng+n);

                obj.(strcat('ind_', idx)) = [obj.(strcat('ind_', idx)), new_indices];
            end
        end
    end
end
