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
        end

        function [] = addVariable(obj, symbolic, type, lb, ub, initial, varargin)
            import casadi.*

            p = inputParser();
            p.FunctionName = 'addVariable';
            
            % TODO: add checks.
            addRequired(p, 'obj');
            addRequired(p, 'symbolic');
            addRequired(p, 'lb');
            addRequired(p, 'ub');
            addRequired(p, 'initial');
            addRequired(p, 'type');
            addOptional(p, 'stage', []);
            addOptional(p, 'sys', []);
            parse(p, obj, symbolic, lb, ub, initial, type, varargin{:})

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
                obj.(strcat('ind_', type)) = new_indices;
            end
            
            % if ismember('stage',p.UsingDefaults)
            %     obj.(strcat('ind_', type)) = [obj.(strcat('ind_', type)), new_indices];
            % else
            %     if ~ismember('sys',p.UsingDefaults)
            %         obj.(strcat('ind_', type)){p.Results.stage, p.Results.sys} = new_indices;
            %     else
            %         obj.(strcat('ind_', type)){p.Results.stage} = new_indices;
            %     end
            % end
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
