classdef Pss < handle
    properties(SetAccess=private)
        regions
        times
    end

    methods
        function obj = Pss(regions, times)
        % Construct an initial PSS active set
        %
        % TODO(@anton) verify regions are nonempty and times are ascending
            arguments
                regions(1,:) cell
                times(1,:) double
            end
            if length(regions) ~= length(times) % TODO(@anton) should we assume t_0 = 0?
                nosnoc.error('size_mismatch', ['Active set components do not have matching size: n_regions=' num2str(length(regions)) ', n_times=' num2str(length(times))]);
            end
            obj.regions = regions;
            obj.times = times;
        end

        function n_steps = get_n_steps(obj)
            n_steps = length(obj.regions);
        end
    end
end
