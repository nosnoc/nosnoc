classdef TestPssToStewart < matlab.unittest.TestCase
    methods(Test)
        function test_basic(tc)
            import casadi.*
            model = nosnoc.model.Pss();
            opts = nosnoc.Options();
            model.x = SX.sym('x', 2);
            model.c = model.x;
            model.S = [1,1;
                -1,1;
                1,-1;
                -1,-1];
            model.F = [1,-1,1,-1;
                1,1,-1,-1];
            model.verify_and_backfill(opts);
            dcs = nosnoc.dcs.Stewart(model);
            dcs.generate_variables(opts);
            dcs.generate_equations(opts);

            tc.assertTrue(isfield(dcs.dims, 'n_theta'));
            tc.assertTrue(isfield(dcs.dims, 'n_lambda'));
            tc.assertTrue(isfield(dcs.dims, 'n_mu'));
            tc.verifyEqual(dcs.dims.n_theta, 4);
            tc.verifyEqual(size(dcs.theta), [4,1]);
            tc.verifyEqual(dcs.dims.n_lambda, 4);
            tc.verifyEqual(size(dcs.lambda), [4,1]);
            tc.verifyEqual(dcs.dims.n_mu, 1);
            tc.verifyEqual(size(dcs.mu), [1,1]);

            % Using string comparison for now.
            % In future it would be nice to be able to check equality
            % via a walk down the expression graph but that is future work.
            tc.verifyEqual(str(dcs.f_x), '[(((theta_1_0-theta_1_1)+theta_1_2)-theta_1_3), (((theta_1_0+theta_1_1)-theta_1_2)-theta_1_3)]');
            tc.verifyEqual(str(dcs.g_ind{1}), '[(-(x_0+x_1)), (x_0-x_1), (x_1-x_0), (x_0+x_1)]');
            tc.verifyEqual(str(dcs.g_alg), '[(mu_1-((x_0+x_1)+lambda_1_0)), (((x_0-x_1)-lambda_1_1)+mu_1), (((x_1-x_0)-lambda_1_2)+mu_1), (((x_0+x_1)-lambda_1_3)+mu_1), ((((theta_1_0+theta_1_1)+theta_1_2)+theta_1_3)-1)]');
        end
        
        function test_simplest(tc)
            import casadi.*

            model = tc.simplest_model();
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            dcs = nosnoc.dcs.Stewart(model);
            dcs.generate_variables(opts);

            tc.verifySize(dcs.theta, [2,1]);
            tc.verifySize(dcs.theta_sys, [1,1]);
            tc.verifySize(dcs.lambda, [2,1]);
            tc.verifySize(dcs.lambda_sys, [1,1]);
            tc.verifySize(dcs.mu, [1,1]);
            tc.verifySize(dcs.mu_sys, [1,1]);
            tc.verifySize(dcs.z_all, [5,1]);

            tc.assertTrue(isfield(dcs.dims, "n_theta"));
            tc.assertTrue(isfield(dcs.dims, "n_lambda"));
            tc.assertTrue(isfield(dcs.dims, "n_mu"));
            tc.verifyEqual(dcs.dims.n_theta, 2);
            tc.verifyEqual(dcs.dims.n_lambda, 2);
            tc.verifyEqual(dcs.dims.n_mu, 1);


            dcs.generate_equations(opts);
            tc.verifySize(dcs.f_x, [1,1]);

            tc.assertTrue(isfield(dcs.dims, 'n_theta'));
            tc.assertTrue(isfield(dcs.dims, 'n_lambda'));
            tc.assertTrue(isfield(dcs.dims, 'n_mu'));
            tc.verifyEqual(dcs.dims.n_theta, 2);
            tc.verifyEqual(size(dcs.theta), [2,1]);
            tc.verifyEqual(dcs.dims.n_lambda, 2);
            tc.verifyEqual(size(dcs.lambda), [2,1]);
            tc.verifyEqual(dcs.dims.n_mu, 1);
            tc.verifyEqual(size(dcs.mu), [1,1]);

            % Using string comparison for now.
            % In future it would be nice to be able to check equality
            % via a walk down the expression graph but that is future work.
            tc.verifyEqual(str(dcs.f_x), '((3*theta_1_0)+theta_1_1)');
            tc.verifyEqual(str(dcs.g_ind{1}), '[x1, (-x1)]');
            tc.verifyEqual(str(dcs.g_alg), '[((x1-lambda_1_0)+mu_1), (mu_1-(x1+lambda_1_1)), ((theta_1_0+theta_1_1)-1)]');
        end

        function test_friction_blocks(tc)
            import casadi.*

            model = tc.friction_blocks();
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            dcs = nosnoc.dcs.Stewart(model);
            dcs.generate_variables(opts);

            tc.verifySize(dcs.theta, [6,1]);
            tc.verifySize(dcs.theta_sys, [1,3]);
            tc.verifySize(dcs.lambda, [6,1]);
            tc.verifySize(dcs.lambda_sys, [1,3]);
            tc.verifySize(dcs.mu, [3,1]);
            tc.verifySize(dcs.mu_sys, [1,3]);
            tc.verifySize(dcs.z_all, [15,1]);

            tc.assertTrue(isfield(dcs.dims, 'n_theta'));
            tc.assertTrue(isfield(dcs.dims, 'n_lambda'));
            tc.assertTrue(isfield(dcs.dims, 'n_mu'));
            tc.verifyEqual(dcs.dims.n_theta, 6);
            tc.verifyEqual(dcs.dims.n_lambda, 6);
            tc.verifyEqual(dcs.dims.n_mu, 3);


            dcs.generate_equations(opts);
            tc.verifySize(dcs.f_x, [7,1]);

            % Using string comparison for now.
            % In future it would be nice to be able to check equality
            % via a walk down the expression graph but that is future work.
            tc.verifyEqual(str(dcs.f_x), '@1=(((q2-q1)-q1)-v1), @2=-0.3, @3=0.3, @4=(((q1-q2)+(q3-q2))-v2), @5=(((q2-q3)-v3)+(10*cos((3.14159*t)))), [((v1*theta_1_0)+(v1*theta_1_1)), ((v2*theta_1_0)+(v2*theta_1_1)), ((v3*theta_1_0)+(v3*theta_1_1)), (((@1+@2)*theta_1_0)+((@1+@3)*theta_1_1)), (((@4*theta_1_0)+(@4*theta_1_1))+((@2*theta_2_0)+(@3*theta_2_1))), (((@5*theta_1_0)+(@5*theta_1_1))+((@2*theta_3_0)+(@3*theta_3_1))), (theta_1_0+theta_1_1)]');
            tc.verifyEqual(str(dcs.g_ind{1}), '[(-v1), v1]');
            tc.verifyEqual(str(dcs.g_ind{2}), '[(-v2), v2]');
            tc.verifyEqual(str(dcs.g_ind{3}), '[(-v3), v3]');
            tc.verifyEqual(str(dcs.g_alg), '@1=1, [(mu_1-(v1+lambda_1_0)), ((v1-lambda_1_1)+mu_1), (mu_2-(v2+lambda_2_0)), ((v2-lambda_2_1)+mu_2), (mu_3-(v3+lambda_3_0)), ((v3-lambda_3_1)+mu_3), ((theta_1_0+theta_1_1)-@1), ((theta_2_0+theta_2_1)-@1), ((theta_3_0+theta_3_1)-@1)]');
        end
    end

    methods
        function model = simplest_model(tc)
            import casadi.*
            model = nosnoc.model.Pss();
            
            x1 = SX.sym('x1');
            model.x = x1;
            model.c = [x1];
            model.S = [-1; 1];

            f_11 = 3;
            f_12 = 1;
            model.F = [f_11, f_12];
            model.x0 = -1;
        end

        function model = friction_blocks(tc)
            import casadi.*
            model = nosnoc.model.Pss();
            model.x0 = [-1;1;-1;-1;1;1;0];
            % differential states
            q1 = SX.sym('q1');
            q2 = SX.sym('q2');
            q3 = SX.sym('q3');
            v1 = SX.sym('v1');
            v2 = SX.sym('v2');
            v3 = SX.sym('v3');
            t = SX.sym('t');
            q = [q1;q2;q3];
            v = [v1;v2;v3];
            model.x = [q;v;t];
            c1 = v1;
            c2 = v2;
            c3 = v3;
            S = [1;-1];
            model.S = {S,S,S};
            model.c = {c1,c2,c3};
            F_external = 0;
            F_input = 10;
            f_base = [v1;...
                v2;...
                v3;...
                (-q1)+(q2-q1)-v1;...
                (q1-q2)+(q3-q2)-v2;...
                (q2-q3)-v3+F_input*cos(pi*t);...
                     ]/6;

            f_base = [v1;...
                v2;...
                v3;...
                (-q1)+(q2-q1)-v1;...
                (q1-q2)+(q3-q2)-v2;...
                (q2-q3)-v3+F_external+F_input*(1*0+1*cos(pi*t));...
                1];
            f_11 = f_base+[0;0;0;-0.3;0;0;0];
            f_12 = f_base+[0;0;0;+0.3;0;0;0];
            f_21 = [0;0;0;0;-0.3;0;0];
            f_22 = [0;0;0;0;0.3;0;0];
            f_31 = [0;0;0;0;0;-0.3;0];
            f_32 = [0;0;0;0;0;0.3;0];
            F1 = [f_11 f_12];
            F2 = [f_21 f_22];
            F3 = [f_31 f_32];

            model.F = {F1 F2 F3};
        end
    end
end
