classdef Pss < handle
    properties(SetAccess=private)
        regions % List of sets of regions which are active. If multiple regions are active, e.g., R_1 and R_2, then a sliding mode between R_1 and R_2 is assumed.
        times   % Switch times, i.e., grid points at which the constant active set changes to a new one, provided in regions.
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
