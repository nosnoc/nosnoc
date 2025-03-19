classdef Pds < nosnoc.activeset.Base
    properties(SetAccess=private)
        active_constraints % List of sets of active constraint indices. E.g. [1,2] means c_{1} = 0 and c_{2} = 0 on an interval specified in times.
        times              % Switch times, i.e., grid points at which the constant active set changes to a new one, provided in regions.
        stages             % Switch (stage, fe) pairs.
    end

    methods
        function obj = Pds(active_constraints, params)
        % Construct an initial Pds active set
        %
        % TODO(@anton) verify regions are nonempty and times are ascending
            arguments
                active_constraints(1,:) cell
                params.times(1,:) double = []
                params.stages(1,:) cell = {}
            end
            if ~isempty(params.times) && length(active_constraints) ~= length(params.times) % TODO(@anton) should we assume t_0 = 0?
                nosnoc.error('size_mismatch', ['Active set components do not have matching size: n_regions=' num2str(length(active_constraints)) ', n_times=' num2str(length(params.times))]);
            elseif ~isempty(params.stages) && length(active_constraints) ~= length(params.stages)
                nosnoc.error('size_mismatch', ['Active set components do not have matching size: n_regions=' num2str(length(active_constraints)) ', n_times=' num2str(length(params.stages))]);
            end
            obj.active_constraints = active_constraints;
            obj.times = params.times;
            obj.stages = params.stages;
        end

        function n_steps = get_n_steps(obj)
            n_steps = length(obj.active_constraints);
        end
    end
end
