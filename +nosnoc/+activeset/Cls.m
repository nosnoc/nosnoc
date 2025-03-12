classdef Cls < handle
    properties(SetAccess=private)
        active_constraints % List of sets of active constraint indices.
        impulses           % List of impulse active sets at switch times.
        times              % Switch times.
        stages             % Switch (stage, fe) pairs.
    end

    methods
        function obj = Cls(active_constraints, impulses, params)
        % Construct an initial Cls active set
        %
        % TODO(@anton) verify regions are nonempty and times are ascending
            arguments
                active_constraints(1,:) cell
                impulses(1,:) logical
                params.times(1,:) double = []
                params.stages(1,:) cell = {}
            end
            nosnoc.error('not_implemented', 'CLS active set passing is not yet implemented')
            if ~isempty(params.times) && length(active_constraints) ~= length(params.times) % TODO(@anton) should we assume t_0 = 0?
                nosnoc.error('size_mismatch', ['Active set components do not have matching size: n_regions=' num2str(length(active_constraints)) ', n_times=' num2str(length(params.times))]);
            elseif ~isempty(params.stages) && length(active_constraints) ~= length(params.stages)
                nosnoc.error('size_mismatch', ['Active set components do not have matching size: n_regions=' num2str(length(active_constraints)) ', n_times=' num2str(length(params.stages))]);
            end
            obj.active_constraints = active_constraints;
            obj.impulses = impulses
            obj.times = params.times;
            obj.stages = params.stages;
        end

        function n_steps = get_n_steps(obj)
            n_steps = length(obj.active_constraints);
        end
    end
end
