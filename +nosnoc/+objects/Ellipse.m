classdef Ellipse < nosnoc.objects.Object
% A shape defined by possibly a union of ellipses (or in 3 dimensions ellipsoids).
    properties(GetAccess=public, SetAccess=private)
        N(1,1) double % double: number of ellipses unioned together.

        x % casadi.SX|casadi.MX: Expression for center of ellipse
        x_dot_lift % casadi.SX|casadi.MX: lifting variables for derivative.

        c % casadi.SX|casadi.MX: Expression for center of ellipse.
        xi % casadi.SX|casadi.MX: Expression for angle of rotation of ellipse.
    end

    properties
        A  % double|cell(double): Matrix that defines the ellipse as x'Ax <= 1
        x0(:,1) double % double: Initial state of the ellipse.
        lbx(:,1) double % double: Lower bound of ellipse state.
        ubx(:,1) double % double: Upper bound of ellipse state.
        x_dot % casadi.SX|casadi.MX: Expression for time derivative of state.
    end

    methods
        function obj = Ellipse(A)
            import casadi.*
            obj.num_objects.count = obj.num_objects.count + 1;
            if ~exist('name')
                obj.name = ['Ellipse' num2str(nosnoc.objects.Ball.num_objects.count)]; 
            end
            obj.A = A;
            
            if iscell(obj.A)
                obj.N = length(obj.A);
            else
                obj.N = 1;
                obj.A = {obj.A};
            end

            obj.n_dim = length(obj.A{1});
            if obj.n_dim > 3 || obj.n_dim < 2
                error("Only 2 or 3 dimensional objects are supported")
            end

            if obj.n_dim == 2
                n_x = 3;
            else
                n_x = 6;
            end
            
            obj.x = SX.sym(['x_' obj.name], n_x);
            if obj.n_dim == 2
                obj.c = obj.x(1:2);
                obj.xi = obj.x(3);
            else
                obj.c = obj.x(1:3);
                obj.xi = obj.x(4:6);
            end
            obj.x_dot = SX(n_x, 1);
            obj.x0 = zeros(n_x,1);
            obj.x_dot_lift = SX.sym(['x_dot_' obj.name], n_x);
            obj.lbx = -inf*ones(n_x,1);
            obj.ubx = inf*ones(n_x,1);
        end

        function poly = to_polygon(obj)
            if obj.n_dim == 2
                for ii=1:length(obj.A)
                    A = obj.A{ii};
                    t=linspace(0,2*pi,360);t(end) = [];
                    uv = [cos(t);sin(t)];
                    Ap = A/A(1,1);
                    c = Ap(1,2);
                    d = sqrt(Ap(2,2));
                    e = sqrt(1/A(1,1));
                    f1 = [e;0];
                    f2 = (e/sqrt(d^2 -c^2))*[-c; 1];
                    F = [f1,f2];
                    xy = F*uv;
                    poly_i = polyshape(xy(1,:), xy(2,:));
                    if ii==1
                        poly = poly_i;
                    else
                        poly = union(poly,poly_i);
                    end
                end
            else
                if length(obj.A) ~= 1
                    warning('Viz of unions of 3d Ellipsoids does not work correctly');
                end
                warning('For now we assume that A is diagonal for vizualization')
                for ii=1:length(obj.A)
                    A = obj.A{ii};
                    [X,Y,Z] = ellipsoid(0,0,0,sqrt(1/A(1,1)),sqrt(1/A(2,2)),sqrt(1/A(3,3)));
                    X = reshape(X, prod(size(X)), 1);
                    Y = reshape(Y, prod(size(Y)), 1);
                    Z = reshape(Z, prod(size(Z)), 1);
                    poly = alphaShape(X, Y, Z, Inf);
                end
            end
        end
    end
end
