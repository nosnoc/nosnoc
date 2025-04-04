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
            n_steps = length(obj.I_0);
        end
    end

    methods(Static)
        function obj = from_pss(pss, model)
        % Construct a heaviside active set from a pss active set.
        % In order to calculate we need the S from the model.
        %
        % Note:
        %      For don't cares the realization of the PSS active set is non-unique.
        %      This means that we make the choice of setting this to the sliding mode, i.e. $\lambda_n = \lambda_p = 0$.
            times = pss.times;
            stages = pss.stages;
            S = model.S{1};
            I_0 = {};
            I_1 = {};
            I_free = {};
            % For each phase we need to translate to step functions
            for ii=1:pss.get_n_steps
                region_ii = pss.regions{ii};
                Sp = S(region_ii, :) == 1;  % For the current region which indicators are positive.
                Sn = S(region_ii, :) == -1; % For the current region which indicators are negative.
                Sdc = S(region_ii, :) == 0; % For the current region which indicators are ignored.

                % Set up initial step activiy
                I_0_ii = Sn(1,:); % Negative means step=0
                I_1_ii = Sp(1,:); % Positive means step=1
                I_free_ii = Sdc(1,:); % Dont care we choose to put in "sliding mode".
                for jj=2:length(region_ii)
                    % If a given index is both positive and negative in two different regions they must be
                    % sliding mode so move that indicator to the I_free set.
                    move_to_free = (I_0_ii & Sp(jj,:)) | I_1_ii & Sn(jj,:);
                    I_0_ii = I_0_ii & ~move_to_free;
                    I_1_ii = I_1_ii & ~move_to_free;
                    I_free_ii = I_free_ii | move_to_free;
                end
                I_0{ii} = find(I_0_ii);
                I_1{ii} = find(I_1_ii);
                I_free{ii} = find(I_free_ii);
            end
            obj = nosnoc.activeset.Heaviside(I_0,I_1,I_free,times=times,stages=stages);
        end
    end
end
