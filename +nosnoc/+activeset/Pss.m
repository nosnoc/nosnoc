classdef Pss < handle
    properties(SetAccess=private)
        regions % List of sets of regions which are active.
        times   % Switch times.
        stages  % Switch (stage, fe) pairs.
    end

    methods
        function obj = Pss(regions, params)
        % Construct an initial PSS active set
        %
        % TODO(@anton) verify regions are nonempty and times are ascending
            arguments
                regions(1,:) cell
                params.times(1,:) double = []
                params.stages(1,:) cell = {} 
            end
            if ~isempty(params.times) && length(regions) ~= length(params.times) % TODO(@anton) should we assume t_0 = 0?
                nosnoc.error('size_mismatch', ['Active set components do not have matching size: n_regions=' num2str(length(regions)) ', n_times=' num2str(length(params.times))]);
            elseif ~isempty(params.stages) && length(regions) ~= length(params.stages)
                nosnoc.error('size_mismatch', ['Active set components do not have matching size: n_regions=' num2str(length(regions)) ', n_times=' num2str(length(params.stages))]);
            end
            obj.regions = regions;
            obj.times = params.times;
            obj.stages = params.stages;
        end

        function n_steps = get_n_steps(obj)
            n_steps = length(obj.regions);
        end
    end
end
