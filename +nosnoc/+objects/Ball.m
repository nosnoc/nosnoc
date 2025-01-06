classdef Ball < nosnoc.objects.Object
% A shape defined by ball in 2 or 3 dimensions.
    properties(GetAccess=public, SetAccess=private)
        c % casadi.SX|casadi.MX: Expression for center of ball

        x % casadi.SX|casadi.MX: Expression for state of ball
    end

    properties
        r(1,1) double {mustBePositive} = 1 % double: Radius of ball.
        x0(:,1) double % double: Initial state of the ball.
        lbx(:,1) double % double: Lower bound of ball state.
        ubx(:,1) double % double: Upper bound of ball state.
        f_rhs % casadi.SX|casadi.MX: Expression for time derivative of state.
        f_rhs_lift % casadi.SX|casadi.MX: lifting variables for derivative.
    end

    methods
        function obj = Ball(r, n_dim, name)
            import casadi.*;
            obj.num_objects.count = obj.num_objects.count + 1;
            obj.timestamp = obj.num_objects.count;
            
            if ~exist('name')
                obj.name = ['Ball' num2str(nosnoc.objects.Ball.num_objects.count)]; 
            end
            
            obj.r = r;
            obj.n_dim = n_dim;

            obj.x = SX.sym(['x_' obj.name], n_dim);
            obj.f_rhs_lift = SX.sym(['f_rhs_' obj.name], n_dim);
            obj.f_rhs = SX(n_dim,1);
            obj.c = obj.x;
            obj.x0 = zeros(n_dim,1);
            obj.lbx = -inf*ones(n_dim,1);
            obj.ubx = inf*ones(n_dim,1);
        end

        function poly = to_polygon(obj)
            if obj.n_dim == 2
                t=linspace(0,2*pi,360).';t(end) = [];
                poly = polyshape(obj.r*cos(t), obj.r*sin(t));
            else
                [x1,y1,z1] = sphere(24);
                x1 = x1(:);
                y1 = y1(:);
                z1 = z1(:);
                P = [x1 y1 z1];
                P = unique(P,'rows');
                poly = alphaShape(P, Inf);
            end
        end
    end
end
