classdef Heaviside < nosnoc.activeset.Base
    properties(SetAccess=private)
        I_0    % List of sets of indices, whose corresponding step function is set to 0.
        I_1    % List of sets of indices, whose corresponding step function is set to 1.
        I_free % List of sets of indices, whose corresponding step function is set to [0,1]. (Corresponds to sliding modes)
        times  % % Switch times, i.e., grid points at which the constant active set changes to a new one.
        stages % Switch (stage, fe) pairs.
    end

    methods
        function obj = Heaviside(I_0, I_1, I_free, params)
        % Construct an initial Heaviside active set
        %
        % TODO(@anton) verify regions are nonempty and times are ascending
            arguments
                I_0(1,:) cell
                I_1(1,:) cell
                I_free(1,:) cell
                params.times(1,:) double = []
                params.stages(1,:) cell = {}
            end
            n_steps = length(I_0);
            if any([length(I_1),length(I_free)] ~= n_steps)
                nosnoc.error('size_mismatch',['Active set components do not have matching size: n_I_0=' num2str(length(I_0)) ', n_I_1=' num2str(length(I_1)), ' n_I_free=' num2str(length(I_free))]);
            end
            if ~isempty(params.times) && n_steps ~= length(params.times) % TODO(@anton) should we assume t_0 = 0?
                nosnoc.error('size_mismatch', ['Active set components do not have matching size: n_regions=' num2str(n_steps) ', n_times=' num2str(length(params.times))]);
            elseif ~isempty(params.stages) && n_steps ~= length(params.stages)
                nosnoc.error('size_mismatch', ['Active set components do not have matching size: n_regions=' num2str(n_steps) ', n_times=' num2str(length(params.stages))]);
            elseif isempty(params.stages) && isempty(params.times)
                nosnoc.error('missing_switches', 'Please pass either the switching times via times, or stages via stages.');
            end
            obj.I_0 = I_0;
            obj.I_1 = I_1;
            obj.I_free = I_free;
            obj.times = params.times;
            obj.stages = params.stages;
        end

        function n_steps = get_n_steps(obj)
            n_steps = length(obj.regions);
        end
    end
end
