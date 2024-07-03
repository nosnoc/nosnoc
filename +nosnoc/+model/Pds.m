classdef Pds < nosnoc.model.Base
% A model that represents a Projected dynamical system. Which has the dynamics:
% $$\dot{x} = \operatorname{P}_{\mathcal{T}_{\mathcal{C}(x)}}(f(x))$$
% with $\operatorname{P}_K(v) = \operatorname{arg\, min}_{s\in K} \frac{1}{2}||s-v||^2_2$,
% and the feasible set set $\mathcal{C} = \{x | c(x) \ge 0\}$.
    properties
        f_x_unconstrained % casadi.SX|casadi.MX: Unconstrained system dynamics expression $f(x)$.

        c % casadi.SX|casadi.MX: The gap functions $c(x)$ used in the definition of the feasible set $\mathcal{C}$. 

        E % double: Square matrix $E \in \mathbb{R}^{n_x\times n_x}$ that is used as the weighting matrix for the projection operator.
    end

    methods
        function obj = Pds()
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

            if size(obj.f_x_unconstrained,1) ~= dims.n_x
                error("nosnoc: f_x_unconstrained has incorrect dimension. It must have the same dimension as x.")
            end

            dims.n_c = size(obj.c,1);

            if ~isempty(obj.E)
                if size(obj.E, 1) ~= dims.n_x
                    error("nosnoc: Projection matrix E must be an n_x by n_x matrix where n_x is the number of functions defining the set C.")
                end
                if size(obj.E, 2) ~= dims.n_x
                    error("nosnoc: Projection matrix E must be an n_x by n_x matrix where n_x is the number of functions defining the set C.")
                end
            else
                obj.E = eye(dims.n_x);
            end

            obj.dims = dims;
        end
    end
end
