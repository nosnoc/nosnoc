classdef NosnocSolver < handle
    properties
        model
        settings
        problem

        solver
        solver_initialization
    end

    methods
        function obj = NosnocSolver(model, settings)
            [obj.solver,obj.solver_initialization,obj.model,obj.settings] = create_nlp_nosnoc(model,settings);
            obj.problem = obj.model.problem;
        end

        function set(obj, type, val)
            if strcmp(type, 'w0')
                % check if val is size w0 then set it
            else
                ind = obj.problem.(strcat('ind_', type));
                flat_ind = sort([ind{:}]);

                if iscell(val)
                    % TODO do automatic fallback on smaller dimension cell arrays
                    if ndims(val) == 2 && size(val, 1) == obj.model.dims.N_stages && size(val,2) == 1
                        for ii=1:obj.model.dims.N_stages
                            for v=ind(ii,:,:)
                                if ~isempty(v) && length(v{1}) == length(val{ii})
                                    obj.solver_initialization.w0(v{1}) = val{ii}; 
                                end
                            end
                        end
                    elseif ndims(val) == 2 && size(val, 1) == obj.model.dims.N_stages && size(val, 2) == obj.model.dims.N_fe
                    end
                else
                    if ndims(val) == 1
                        if length(val) == length(ind)
                            obj.solver_initialization.w0(flat_ind) = val;
                        end
                    else
                        error('nosnoc: set should be a cell array or a flat array')
                    end
                end
            end
        end
        
        function varargout = solve(obj)
            
        end 
    end
end

