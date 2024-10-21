classdef TestPssToHeaviside < matlab.unittest.TestCase
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
            dcs = nosnoc.dcs.Heaviside(model);
            dcs.generate_variables(opts);
            dcs.generate_equations(opts);

            tc.assertTrue(isfield(dcs.dims, 'n_alpha'));
            tc.assertTrue(isfield(dcs.dims, 'n_lambda'));
            tc.verifyEqual(dcs.dims.n_alpha, 2);
            tc.verifyEqual(size(dcs.alpha), [2,1]);
            tc.verifyEqual(dcs.dims.n_lambda, 2);
            tc.verifyEqual(size(dcs.lambda_p), [2,1]);
            tc.verifyEqual(size(dcs.lambda_n), [2,1]);

            % Using string comparison for now.
            % In future it would be nice to be able to check equality
            % via a walk down the expression graph but that is future work.
            tc.verifyEqual(str(dcs.f_x), '@1=(alpha_1_0*alpha_1_1), @2=1, @3=((@2-alpha_1_0)*alpha_1_1), @4=(alpha_1_0*(@2-alpha_1_1)), @5=((@2-alpha_1_0)*(@2-alpha_1_1)), [(((@1-@3)+@4)-@5), (((@1+@3)-@4)-@5)]');
            tc.verifyEqual(str(dcs.g_alg), '[((x_0-lambda_p_1_0)+lambda_n_1_0), ((x_1-lambda_p_1_1)+lambda_n_1_1)]');
        end
        
        function test_simplest(tc)
            import casadi.*

            model = tc.simplest_model();
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            dcs = nosnoc.dcs.Heaviside(model);
            dcs.generate_variables(opts);

            tc.assertTrue(isfield(dcs.dims, 'n_alpha'));
            tc.assertTrue(isfield(dcs.dims, 'n_lambda'));
            tc.verifyEqual(dcs.dims.n_alpha, 1);
            tc.verifyEqual(dcs.dims.n_lambda, 1);

            tc.verifySize(dcs.alpha, [1,1]);
            tc.verifySize(dcs.lambda_n, [1,1]);
            tc.verifySize(dcs.lambda_p, [1,1]);


            dcs.generate_equations(opts);
            tc.verifySize(dcs.f_x, [1,1]);


            % Using string comparison for now.
            % In future it would be nice to be able to check equality
            % via a walk down the expression graph but that is future work.
            tc.verifyEqual(str(dcs.f_x), '((3*(1-alpha_1))+alpha_1)');
            tc.verifyEqual(str(dcs.g_alg), '((x1-lambda_p_1)+lambda_n_1)');
        end

        function test_friction_blocks(tc)
            import casadi.*

            model = tc.friction_blocks();
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            dcs = nosnoc.dcs.Heaviside(model);
            dcs.generate_variables(opts);

            tc.verifySize(dcs.alpha, [3,1]);
            tc.verifySize(dcs.lambda_p, [3,1]);
            tc.verifySize(dcs.lambda_n, [3,1]);

            tc.assertTrue(isfield(dcs.dims, 'n_alpha'));
            tc.assertTrue(isfield(dcs.dims, 'n_lambda'));
            tc.verifyEqual(dcs.dims.n_alpha, 3);
            tc.verifyEqual(dcs.dims.n_lambda, 3);


            dcs.generate_equations(opts);
            tc.verifySize(dcs.f_x, [7,1]);

            % Using string comparison for now.
            % In future it would be nice to be able to check equality
            % via a walk down the expression graph but that is future work.
            tc.verifyEqual(str(dcs.f_x), '@1=1, @2=(@1-alpha_1), @3=(((q2-q1)-q1)-v1), @4=-0.3, @5=0.3, @6=(((q1-q2)+(q3-q2))-v2), @7=(((q2-q3)-v3)+(10*cos((3.14159*t)))), [((v1*alpha_1)+(v1*@2)), ((v2*alpha_1)+(v2*@2)), ((v3*alpha_1)+(v3*@2)), (((@3+@4)*alpha_1)+((@3+@5)*@2)), (((@6*alpha_1)+(@6*@2))+((@4*alpha_2)+(@5*(@1-alpha_2)))), (((@7*alpha_1)+(@7*@2))+((@4*alpha_3)+(@5*(@1-alpha_3)))), @1]');
            tc.verifyEqual(str(dcs.g_alg), '[((v1-lambda_p_1)+lambda_n_1), ((v2-lambda_p_2)+lambda_n_2), ((v3-lambda_p_3)+lambda_n_3)]');
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
