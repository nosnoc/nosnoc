classdef TestPssStewart < matlab.unittest.TestCase
    methods(Test)
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
            % TODO probably verify some other things about f_x et al.
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
