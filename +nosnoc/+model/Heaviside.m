classdef Heaviside < nosnoc.model.Base
% A nonsmooth model which allows a more general than a PSS as it implements the Aizermann-Pyatnitskii extension.
% In particular it allows multipliers for example in the form :math:`1-\alpha_1\alpha_2`.
    properties
        f_x % casadi.SX|casadi.MX: $\xdot = f(x,\alpha)$ i.e. the state dynamics. 

        c % casadi.SX|casadi.MX: $c(x)$ representing the switching function for the $\alpha$ multipliers. 

        % casadi.SX|casadi.MX: Step functions $\alpha$ which are definied by: $\alpha = 1$ if $c(x) > 0$,
        % $\alpha = 0$ if $c(x) < 0$, and $\alpha \in [0, 1]$ if $c(x) = 0$.
        alpha
    end

    methods
        function obj = Heaviside()
            obj = obj@nosnoc.model.Base();
        end
        
        function verify_and_backfill(obj, opts)
            arguments
                obj
                opts nosnoc.Options
            end
            import casadi.*
            verify_and_backfill@nosnoc.model.Base(obj,opts);

            dims = obj.dims;

            % check how many subsystems are present
            if iscell(obj.c)
                dims.n_sys = length(obj.c);
            else
                obj.c = {obj.c};
                dims.n_sys = 1;
            end

            dims.n_c_sys = [];
            for ii = 1:dims.n_sys
                % dimensions of c
                dims.n_c_sys  = [dims.n_c_sys;length(obj.c{ii})];
            end

            if size(obj.alpha,1) ~= sum(dims.n_c_sys)
                error('nosnoc: There needs to be a step function alpha for each switching function in c');
            end
            dims.n_alpha = size(obj.alpha, 1);

            obj.dims = dims;
        end
    end
end
