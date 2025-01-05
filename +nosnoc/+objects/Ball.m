classdef Ball < nosnoc.objects.Object
% A shape defined by possibly a union of half-space convex polygons.
    properties(GetAccess=public, SetAccess=private)
        c

        x
    end

    properties
        r
        x0
        lbx
        ubx
        x_dot
        x_dot_lift
    end

    methods
        function obj = Ball(r, n_dim, name)
            import casadi.*;
            obj.num_objects.count = obj.num_objects.count + 1;
            if ~exist('name')
                obj.name = ['Ball' num2str(nosnoc.objects.Ball.num_objects.count)]; 
            end
            
            obj.r = r;
            obj.n_dim = n_dim;

            obj.x = SX.sym(['x_' obj.name], n_dim);
            obj.x_dot_lift = SX.sym(['x_dot_' obj.name], n_dim);
            obj.x_dot = SX(n_dim,1);
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
