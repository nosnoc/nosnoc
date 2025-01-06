% BSD 2-Clause License

% Copyright (c) 2023, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% This file is part of NOSNOC.
classdef PDSObjects < nosnoc.model.Base
    properties (Access=public)
        
        % Algebraics
        lambda % casadi.SX|casadi.MX: Standard PDS lagrange multipliers.
        mu % casadi.SX|casadi.MX: Lagrange multipliers for SDF KKT conditions.
        alpha % casadi.SX|casadi.MX: SDF objectives.
        p_d % casadi.SX|casadi.MX: Contact point coordinates.
        y1_d % casadi.SX|casadi.MX: Auxiliary points for padded polytopes.
        y2_d % casadi.SX|casadi.MX: Auxiliary points for padded polytopes.
        lambda_t % casadi.SX|casadi.MX: Multiplier for stick slip friction.
        v_t % casadi.SX|casadi.MX: Lifted tangential velocity.
        normal_lift % casadi.SX|casadi.MX: Lifted normal.
        gamma_f % casadi.SX|casadi.MX: Friction slack variable

        f_rhs_lift % casadi.SX|casadi.MX: Lifted f_rhs

        G_friction % casadi.SX|casadi.MX: Expression for complementarity values G for friction. 
        H_friction % casadi.SX|casadi.MX: Expression for complementarity values H for friction.
        g_friction % casadi.SX|casadi.MX: Expression for equality constraints for friction

        % r = 0.01 % Padding constant
        
        g_d % casadi.SX|casadi.MX: Constraint expressions for SDF
        g_kkt % casadi.SX|casadi.MX: KKT expression for SDF
        normal % casadi.SX|casadi.MX: expression for normal vector to contact.
        normal0 % double: initial values of normal vector
        tangent % casadi.SX|casadi.MX: expression for tangent vector to contact
        tangent0 % double: initial values of tangent vector

        % Distance functions
        c % casadi.SX|casadi.MX: Standard PDS distance functions.

        % objects
        objects % nosnoc.objects.Object: List of objects.
        contacts % cell: List of contact pairs.

    end

    methods (Access=public)
        function obj = PDSObjects()
            dims = struct;
        end

        function verify_and_backfill(obj, opts)
            import casadi.*
            verify_and_backfill@nosnoc.model.Base(obj,opts);

            dims = obj.dims;
            
            dims.n_c = size(obj.c, 1);
            dims.n_lambda = dims.n_c;
            dims.n_y1d = size(obj.y1_d, 1);
            dims.n_y2d = size(obj.y2_d, 1);
            dims.n_mu = size(obj.mu, 1);
            dims.n_pd = size(obj.p_d,1);
            dims.n_alpha = size(obj.alpha,1);
            dims.n_tangent = size(obj.tangent,1);
            dims.n_normal = size(obj.normal,1);
            dims.n_normal_lift = size(obj.normal_lift,1);
            dims.n_v_t = size(obj.v_t,1);
            dims.n_gamma_f = size(obj.gamma_f,1);
            dims.n_lambda_t = size(obj.lambda_t,1);
            obj.dims = dims;

        end

        function addContact(obj, shape1, shape2, mu)
        % Add contact between two shapes of compatible dimension possibly with "friction" model.
            arguments
                obj nosnoc.model.PDSObjects
                shape1 nosnoc.objects.Object % First shape
                shape2 nosnoc.objects.Object % Second shape
                mu = 0 % Friction coef
            end
            if shape1.n_dim ~= shape2.n_dim
                error('Dimensions mismatch');
            end
            
            if isa(shape1, 'nosnoc.objects.Ball') && isa(shape2, 'nosnoc.objects.Ball')
                obj.addBallBall(shape1, shape2, mu);
            elseif isa(shape1, 'nosnoc.objects.Ellipse') && isa(shape2, 'nosnoc.objects.Ellipse')
                obj.addEllipseEllipse(shape1, shape2, mu);
            elseif isa(shape1, 'nosnoc.objects.Ball') && isa(shape2, 'nosnoc.objects.Ellipse')
                obj.addEllipseBall(shape2,shape1, mu);
            elseif isa(shape1, 'nosnoc.objects.Ellipse') && isa(shape2, 'nosnoc.objects.Ball')
                obj.addEllipseBall(shape1,shape2, mu);
            else
                error('Unknown shape types');
            end

            if ~ismember(shape1, obj.objects)
                obj.objects = [obj.objects, shape1];
            end

            if ~ismember(shape2, obj.objects)
                obj.objects = [obj.objects, shape2];
            end

            [~,ind] = sort([obj.objects.timestamp]);
            obj.objects = obj.objects(ind);

            obj.contacts{end+1} = {shape1, shape2, mu};

            obj.x = [];
            obj.x0 = [];
            obj.lbx = [];
            obj.ubx = [];
            for ii=1:length(obj.objects)
                obj.x = vertcat(obj.x, obj.objects(ii).x);
                obj.x0 = vertcat(obj.x0, obj.objects(ii).x0);
                obj.lbx = vertcat(obj.lbx, obj.objects(ii).lbx);
                obj.ubx = vertcat(obj.ubx, obj.objects(ii).ubx);
            end
        end
    end

    methods (Access=private)
        function addBallBall(obj, ball1, ball2, mu_f)
        % Function adds a contact between two balls in 2 or 3 dimensions.
        % This can be done without an implcit distance function and so we do it directly 
        % as an explicit function in $x$.
            import casadi.*

            % Lambda and gap functions
            lambda = SX.sym(['lambda_' ball1.name '_' ball2.name]);
            c = sum((ball1.x - ball2.x).^2) - (ball1.r + ball2.r)^2;

            % Get normals
            nabla_c_x1 = c.jacobian(ball1.x)';
            nabla_c_x2 = c.jacobian(ball2.x)';

            % Update dynamics of both balls with normal times lagrange multiplier
            ball1.f_rhs = ball1.f_rhs + nabla_c_x1*lambda;
            ball2.f_rhs = ball2.f_rhs + nabla_c_x2*lambda;

            % Normal expression is just nabla_c as c is explicit in this case.
            normal = [nabla_c_x1;nabla_c_x2];
            obj.normal = [obj.normal;normal];
            obj.normal0 = [obj.normal0;0;1;0;-1];

            % Handle friction if necessary.
            if mu_f
                if ball1.n_dim == 2
                    % Create tangent variables
                    tangent = SX.sym('tangent', 2);
                    obj.tangent = [obj.tangent; tangent];
                    obj.tangent0 = [obj.tangent0;1;0];
                    % Create tangential velocity variables.
                    v_t = SX.sym('v_t', 2);
                    obj.v_t = [obj.v_t; v_t];
                    % Friction values
                    lambda_t = SX.sym('lambda_t', 2);
                    obj.lambda_t = [obj.lambda_t; lambda_t];
                    % Friction slacks
                    gamma_f = SX.sym('gamma_f', 1);
                    obj.gamma_f = [obj.gamma_f; gamma_f];
                    % lifted normal variables.
                    normal_lift = SX.sym('normal_lift', 4);
                    obj.normal_lift = [obj.normal_lift; normal_lift];

                    ball1.f_rhs = ball1.f_rhs + lambda_t(1)*tangent - lambda_t(2)*tangent;
                    ball2.f_rhs = ball2.f_rhs - lambda_t(1)*tangent + lambda_t(2)*tangent;
                    
                    g_friction = [normal_lift - normal; % Lift normal (helps with convergence) TODO(@anton) allow for not lifting
                        tangent - [-normal_lift(2);normal_lift(1)]; % get tangent (lifted again) TODO(@anton) allow for not lifting
                        v_t(1) - (dot(ball1.f_rhs_lift, tangent)/dot(tangent,tangent)); % Tangent velocities of both objects 
                        v_t(2) - (dot(ball2.f_rhs_lift, tangent)/dot(tangent,tangent))];
                    obj.g_friction = [obj.g_friction; g_friction];

                    % Complementarity constraints for friction in the positive or negative tangent direction must be zero
                    % also slack that relates the normal and tangential "forces"
                    G_friction = [v_t(1) - v_t(2) + gamma_f;
                        v_t(2) - v_t(1) + gamma_f;
                        mu_f*lambda - lambda_t(1) - lambda_t(2);
                        lambda_t(1)]; % TODO(@anton) can we avoid this as it violates MPCC-LICQ
                    obj.G_friction = [obj.G_friction; G_friction];
                    H_friction = [lambda_t; gamma_f;lambda_t(2)];
                    obj.H_friction = [obj.H_friction; H_friction];
                    
                else
                    error("3d friction not yet supported")
                end
            end
            
            % update model
            obj.c = vertcat(obj.c,c);
            obj.lambda = vertcat(obj.lambda, lambda);
        end

        function addEllipseBall(obj, ellipse, ball, mu_f)
        % Function adds a contact between a circle and an ellipse in 2d or a ball and
        % an ellipsoid in 3 dimensions.
        % This is done through implicit signed distance functions. See :cite:`Tracy2023`
        % and :cite:`Dietz2024` for more details.
            import casadi.*
            n_dim = ball.n_dim;

            % Helper Functions
            if n_dim == 2
                % Generate 2d rotation matrix function.
                theta = SX.sym('theta');
                Rot = [cos(theta) -sin(theta);...
                    sin(theta) cos(theta)];
                R = Function('R', {theta}, {Rot});
                % Generate vector cross product function
                a = SX.sym('a',2);
                b = SX.sym('b',2);
                cross_fun = Function('cross', {a,b}, {a(1)*b(2) - a(2)*b(1)});
            else
                % Generate 3d rotation matrix function.
                rx = SX.sym('rx',1);
                ry = SX.sym('ry',1);
                rz = SX.sym('rz',1);
                Rot = [cos(ry)*cos(rz), sin(rx)*sin(ry)*cos(rz) - cos(rx)*sin(rz),  cos(rx)*sin(ry)*cos(rz) + sin(rx)*sin(rz);
                    cos(ry)*sin(rz), sin(rx)*sin(ry)*sin(rz) + cos(rx)*cos(rz),  cos(rx)*sin(ry)*sin(rz) - sin(rx)*cos(rz);
                    -sin(ry), sin(rx)*cos(ry),  cos(rx)*cos(ry)];
                R = Function('R', {[rx;ry;rz]}, {Rot});

                % Generate vector cross product function
                a = SX.sym('a',3);
                b = SX.sym('b',3);
                cross_fun = Function('cross', {a,b}, {cross(a,b)});
            end

            % Define contacts for each subellipse in the ellipse object
            for ii=1:ellipse.N
                % Generate required variables.
                lambda = SX.sym(['lambda_' ellipse.name '_' ball.name]);
                obj.lambda = vertcat(obj.lambda, lambda);
                alpha = SX.sym(['alpha_' ellipse.name '_' ball.name]);
                obj.alpha = vertcat(obj.alpha, alpha);
                p_d = SX.sym(['p_d_' ellipse.name '_' ball.name], n_dim);
                obj.p_d = vertcat(obj.p_d, p_d);
                mu = SX.sym(['mu_' ellipse.name '_' ball.name], 2);
                obj.mu = vertcat(obj.mu, mu);

                % Get ellipse matrix and ball radius
                A = ellipse.A{ii};
                r = ball.r;

                % Get KKT multipliers for each object
                mu_1c = mu(end-1);
                mu_2c = mu(end);

                % expanding constraints for ellipse and ball.
                g_d = [(p_d-ellipse.c)'*R(ellipse.xi)*A*R(ellipse.xi)'*(p_d-ellipse.c) - alpha;
                      (1/r^2)*(p_d-ball.c)'*(p_d-ball.c) - alpha];
                obj.g_d = vertcat(obj.g_d, g_d);

                % Lagrangian for the implicit SDF
                L = alpha - 1 + mu'*g_d;

                % KKT stationarity functions for SDF
                % g_kkt = [1 - mu_1c - mu_2c;
                %     mu_1c*(R(ellipse.xi)*(A + A')*R(ellipse.xi)')*(p_d-ellipse.c) + mu_2c*(2/r^2)*(p_d-ball.c)];
                g_kkt = L.jacobian([alpha;p_d])';
                obj.g_kkt = vertcat(obj.g_kkt, g_kkt);

                % normal expression (nabla g_d with respect to positions)
                ntr = 2/ball.r^2*(p_d - ball.c)*mu_2c;
                normal_ball = [-ntr];
                normal_ellipse =[ntr;
                    cross_fun(p_d-ellipse.c, ntr)];

                obj.normal = [obj.normal;normal_ball;normal_ellipse];
                obj.normal0 = [obj.normal0;0;1;0;-1;0];

                % Update dynamics of both objects
                ellipse.f_rhs = ellipse.f_rhs + normal_ellipse*lambda;
                ball.f_rhs = ball.f_rhs + normal_ball*lambda;

                % Do friction. See above for explanation.
                if mu_f
                    if ball.n_dim == 2 
                        tangent = SX.sym('tangent', 2);
                        obj.tangent = [obj.tangent; tangent];
                        obj.tangent0 = [obj.tangent0;1;0];
                        v_t = SX.sym('v_t', 2);
                        obj.v_t = [obj.v_t; v_t];
                        lambda_t = SX.sym('lambda_t', 2);
                        obj.lambda_t = [obj.lambda_t; lambda_t];
                        gamma_f = SX.sym('gamma_f', 1);
                        obj.gamma_f = [obj.gamma_f; gamma_f];
                        normal_lift = SX.sym('normal_lift', 4);
                        obj.normal_lift = [obj.normal_lift; normal_lift];
                        
                        ellipse.f_rhs = ellipse.f_rhs + lambda_t(1)*[tangent;cross_fun(p_d-ellipse.c, tangent)] - lambda_t(2)*[tangent;cross_fun(p_d-ellipse.c, tangent)];
                        ball.f_rhs = ball.f_rhs - lambda_t(1)*tangent + lambda_t(2)*tangent;
                        
                        g_friction = [normal_lift - [normal_ball;normal_ellipse(1:2)]; % Lift normal (helps with convergence) TODO(@anton) allow for not lifting

                            tangent - [-normal_lift(2);normal_lift(1)]; % get tangent (lifted again) TODO(@anton) allow for not lifting
                            v_t(1) - (dot(ellipse.f_rhs_lift(1:2), tangent)/dot(tangent,tangent) + ellipse.f_rhs_lift(3)*sqrt(max(sum((p_d - ellipse.c).^2),1e-9))); % Tangent velocities of both objects 
                            v_t(2) - (dot(ball.f_rhs_lift, tangent)/dot(tangent,tangent))];
                        obj.g_friction = [obj.g_friction; g_friction];
                        
                        G_friction = [v_t(1) - v_t(2) + gamma_f;
                            v_t(2) - v_t(1) + gamma_f;
                            mu_f*lambda - lambda_t(1) - lambda_t(2)]; % TODO(@anton) can we avoid this as it violates MPCC-LICQ
                        obj.G_friction = [obj.G_friction; G_friction];
                        H_friction = [lambda_t; gamma_f];
                        obj.H_friction = [obj.H_friction; H_friction];
                    else
                        error("3d friction not yet supported")
                    end
                end
                
                obj.c = vertcat(obj.c, alpha - 1);
            end
        end
        
        function addEllipseEllipse(obj, ellipse1, ellipse2, mu_f)
        % Function adds a contact between two ellipses in 2d or two ellipsoids in 3d.
        % This is done through implicit signed distance functions. See :cite:`Tracy2023`
        % and :cite:`Dietz2024` for more details.
            import casadi.*
            n_dim = ellipse2.n_dim;

            % Helper Functions
            if n_dim == 2
                % Generate 2d rotation matrix function.
                theta = SX.sym('theta');
                Rot = [cos(theta) -sin(theta);...
                    sin(theta) cos(theta)];
                R = Function('R', {theta}, {Rot});
                
                % Generate vector cross product function
                a = SX.sym('a',2);
                b = SX.sym('b',2);
                cross_fun = Function('cross', {a,b}, {a(1)*b(2) - a(2)*b(1)});
            else
                % Generate 3d rotation matrix function.
                rx = SX.sym('rx',1);
                ry = SX.sym('ry',1);
                rz = SX.sym('rz',1);
                Rot = [cos(ry)*cos(rz), sin(rx)*sin(ry)*cos(rz) - cos(rx)*sin(rz),  cos(rx)*sin(ry)*cos(rz) + sin(rx)*sin(rz);
                    cos(ry)*sin(rz), sin(rx)*sin(ry)*sin(rz) + cos(rx)*cos(rz),  cos(rx)*sin(ry)*sin(rz) - sin(rx)*cos(rz);
                    -sin(ry), sin(rx)*cos(ry),  cos(rx)*cos(ry)];
                R = Function('R', {[rx;ry;rz]}, {Rot});
                
                % Generate vector cross product function
                a = SX.sym('a',3);
                b = SX.sym('b',3);
                cross_fun = Function('cross', {a,b}, {cross(a,b)});
            end

            % Define contacts for each pair of subellipses
            for ii=1:ellipse1.N
                for jj=1:ellipse2.N
                    % Generate required variables for 
                    lambda = SX.sym(['lambda_' ellipse1.name '_' ellipse2.name]);
                    obj.lambda = vertcat(obj.lambda, lambda);
                    alpha = SX.sym(['alpha_' ellipse1.name '_' ellipse2.name]);
                    obj.alpha = vertcat(obj.alpha, alpha);
                    p_d = SX.sym(['p_d_' ellipse1.name '_' ellipse2.name], n_dim);
                    obj.p_d = vertcat(obj.p_d, p_d);
                    mu = SX.sym(['mu_' ellipse1.name '_' ellipse2.name], 2);
                    obj.mu = vertcat(obj.mu, mu);

                    % Get corresponding ellipse matrices
                    A1 = ellipse1.A{ii};
                    A2 = ellipse2.A{jj};

                    % Get KKT multipliers for each object
                    mu_1c = mu(end-1);
                    mu_2c = mu(end);

                    % expanding constraints for ellipses defining $\mathcal{S}(x)$.
                    g_d = [(p_d-ellipse1.c)'*R(ellipse1.xi)*A1*R(ellipse1.xi)'*(p_d-ellipse1.c) - alpha;
                        (p_d-ellipse2.c)'*R(ellipse2.xi)*A2*R(ellipse2.xi)'*(p_d-ellipse2.c) - alpha];
                    obj.g_d = vertcat(obj.g_d, g_d);

                    % Lagrangian for the implicit SDF
                    L = alpha - 1 + mu'*g_d;

                    % KKT stationarity functions for SDF
                    % g_kkt = [1 - mu_1c - mu_2c;
                    %     mu_1c*(R(ellipse1.xi)*(A + A')*R(ellipse1.xi)')*(p_d-ellipse1.c) + mu_2c*(2/r^2)*(p_d-ellipse2.c)];
                    g_kkt = L.jacobian([alpha;p_d])';
                    obj.g_kkt = vertcat(obj.g_kkt, g_kkt);

                    % normal expression (nabla g_d with respect to positions)
                    ntr = -2*R(ellipse1.xi)*A1*R(ellipse1.xi)'*(p_d - ellipse1.c)*mu_1c;
                    normal_ellipse1 =[ntr;
                        cross_fun(p_d-ellipse1.c, ntr)];
                    normal_ellipse2 = [-ntr;
                        -cross_fun(p_d-ellipse2.c, ntr)];

                    obj.normal = [obj.normal;normal_ellipse1;normal_ellipse2];
                    obj.normal0 = [obj.normal0;0;1;0;-1;0];

                    % Update dynamics of both objects
                    ellipse1.f_rhs = ellipse1.f_rhs + normal_ellipse1*lambda;
                    ellipse2.f_rhs = ellipse2.f_rhs + normal_ellipse2*lambda;

                    % Do friction. See above for explanation.
                    if mu_f
                        if ellipse2.n_dim == 2 %TODO does this need changes?
                            tangent = SX.sym('tangent', 2);
                            obj.tangent = [obj.tangent; tangent];
                            obj.tangent0 = [obj.tangent0;1;0];
                            v_t = SX.sym('v_t', 2);
                            obj.v_t = [obj.v_t; v_t];
                            lambda_t = SX.sym('lambda_t', 2);
                            obj.lambda_t = [obj.lambda_t; lambda_t];
                            gamma_f = SX.sym('gamma_f', 1);
                            obj.gamma_f = [obj.gamma_f; gamma_f];
                            normal_lift = SX.sym('normal_lift', 4);
                            obj.normal_lift = [obj.normal_lift; normal_lift];
                            
                            ellipse1.f_rhs = ellipse1.f_rhs + lambda_t(1)*[tangent;cross_fun(p_d-ellipse1.c, tangent)] - lambda_t(2)*[tangent;cross_fun(p_d-ellipse1.c, tangent)];
                            ellipse2.f_rhs = ellipse2.f_rhs - lambda_t(1)*tangent + lambda_t(2)*tangent;
                            
                            g_friction = [normal_lift - [normal_ellipse2;normal_ellipse1(1:2)]; % Lift normal (helps with convergence) TODO(@anton) allow for not lifting
                                tangent - [-normal_lift(2);normal_lift(1)]; % get tangent (lifted again) TODO(@anton) allow for not lifting
                                v_t(1) - (dot(ellipse1.f_rhs_lift(1:2), tangent)/dot(tangent,tangent) + ellipse1.f_rhs_lift(3)*norm(p_d - ellipse1.c)); % Tangent velocities of both objects 
                                v_t(2) - (dot(ellipse2.f_rhs_lift, tangent)/dot(tangent,tangent))];
                            obj.g_friction = [obj.g_friction; g_friction];
                            
                            G_friction = [v_t(1) - v_t(2) + gamma_f;
                                v_t(2) - v_t(1) + gamma_f;
                                mu_f*lambda - lambda_t(1) - lambda_t(2);
                                lambda_t(1);
                                lambda]; % TODO(@anton) can we avoid this as it violates MPCC-LICQ
                            obj.G_friction = [obj.G_friction; G_friction];
                            H_friction = [lambda_t; gamma_f;lambda_t(2);gamma_f];
                            obj.H_friction = [obj.H_friction; H_friction];
                        else
                            error("3d friction not yet supported")
                        end
                    end
                    
                    obj.c = vertcat(obj.c, alpha - 1);
                end
            end
        end
    end
end
