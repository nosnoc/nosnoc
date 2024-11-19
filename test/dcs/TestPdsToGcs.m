classdef TestPdsToGcs < matlab.unittest.TestCase
    methods(Test)
        function test_basic(tc)
            import casadi.*

            model = tc.simplest_model();
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            dcs = nosnoc.dcs.Gcs(model);
            dcs.generate_variables(opts);

            tc.verifySize(dcs.lambda, [1,1]);
            tc.verifySize(dcs.c_lift, [1,1]);

            tc.assertTrue(isfield(dcs.dims, "n_c"));
            tc.assertTrue(isfield(dcs.dims, "n_c_lift"));
            tc.assertTrue(isfield(dcs.dims, "n_lambda"));
            tc.verifyEqual(dcs.dims.n_c, 1);
            tc.verifyEqual(dcs.dims.n_c_lift, 1);
            tc.verifyEqual(dcs.dims.n_lambda, 1);
            
            dcs.generate_equations(opts);
            % Using string comparison for now.
            % In future it would be nice to be able to check equality
            % via a walk down the expression graph but that is future work.
            tc.verifyEqual(str(dcs.f_x), '[x_1, (lambda-x_0)]');
            tc.verifyEqual(str(dcs.g_c_lift), '(c_lift-(x_1+0.2))');
            tc.verifyEqual(str(dcs.c), 'c_lift');
        end

        function test_ocp(tc)
            import casadi.*

            model = tc.pds_ocp();
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            dcs = nosnoc.dcs.Gcs(model);
            dcs.generate_variables(opts);

            tc.verifySize(dcs.lambda, [1,1]);
            tc.verifySize(dcs.c_lift, [1,1]);

            tc.assertTrue(isfield(dcs.dims, "n_c"));
            tc.assertTrue(isfield(dcs.dims, "n_c_lift"));
            tc.assertTrue(isfield(dcs.dims, "n_lambda"));
            tc.verifyEqual(dcs.dims.n_c, 1);
            tc.verifyEqual(dcs.dims.n_c_lift, 1);
            tc.verifyEqual(dcs.dims.n_lambda, 1);
            
            dcs.generate_equations(opts);
            % Using string comparison for now.
            % In future it would be nice to be able to check equality
            % via a walk down the expression graph but that is future work.
            tc.verifyEqual(str(dcs.f_x), '@1=0.25, @2=0.0625, [(((-0.2*sq((x_0+1)))+u1_0)+(((@1*x_0)+(@1*x_0))*lambda)), (((-0.4*(x_1+3))+u1_1)+(((@2*x_1)+(@2*x_1))*lambda))]');
            tc.verifyEqual(str(dcs.g_c_lift), '(c_lift-((((0.25*x_0)*x_0)+((0.0625*x_1)*x_1))-1))');
            tc.verifyEqual(str(dcs.c), 'c_lift');
        end

        function test_open_gate(tc)
            import casadi.*

            model = tc.open_gate();
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            dcs = nosnoc.dcs.Gcs(model);
            dcs.generate_variables(opts);

            tc.verifySize(dcs.lambda, [4,1]);
            tc.verifySize(dcs.c_lift, [4,1]);

            tc.assertTrue(isfield(dcs.dims, "n_c"));
            tc.assertTrue(isfield(dcs.dims, "n_c_lift"));
            tc.assertTrue(isfield(dcs.dims, "n_lambda"));
            tc.verifyEqual(dcs.dims.n_c, 4);
            tc.verifyEqual(dcs.dims.n_c_lift, 4);
            tc.verifyEqual(dcs.dims.n_lambda, 4);
            
            dcs.generate_equations(opts);
            % Using string comparison for now.
            % In future it would be nice to be able to check equality
            % via a walk down the expression graph but that is future work.
            tc.verifyEqual(str(dcs.f_x), '@1=(x3_0-x1_0), @2=(x3_1-x1_1), @3=sqrt((sq(@1)+sq(@2))), @4=(x2_0-x1_0), @5=(x2_1-x1_1), @6=sqrt((sq(@4)+sq(@5))), @7=(x3_0-x2_0), @8=(x3_1-x2_1), @9=sqrt((sq(@7)+sq(@8))), [(u1_0-(((@1/@3)*lambda_0)+((@4/@6)*lambda_2))), (u1_1+(lambda_3-(((@2/@3)*lambda_0)+((@5/@6)*lambda_2)))), (u2_0+(((@4/@6)*lambda_2)-((@7/@9)*lambda_1))), (u2_1+(((@5/@6)*lambda_2)-((@8/@9)*lambda_1))), (((@1/@3)*lambda_0)+((@7/@9)*lambda_1)), (((@2/@3)*lambda_0)+((@8/@9)*lambda_1)), (-lambda_3)]');
            tc.verifyEqual(str(dcs.g_c_lift), '@1=2, [(c_lift_0-(sqrt((sq((x3_0-x1_0))+sq((x3_1-x1_1))))-@1)), (c_lift_1-(sqrt((sq((x3_0-x2_0))+sq((x3_1-x2_1))))-@1)), (c_lift_2-(sqrt((sq((x2_0-x1_0))+sq((x2_1-x1_1))))-@1)), (c_lift_3-((x1_1-gate)-1))]');
            tc.verifyEqual(str(dcs.c), '[c_lift_0, c_lift_1, c_lift_2, c_lift_3]');
        end
    end

    methods
        function model = simplest_model(tc)
            import casadi.*
            model = nosnoc.model.Pds();

            x = SX.sym('x',2);
            model = nosnoc.model.Pds();
            model.x0 = [sqrt(2)/2;sqrt(2)/2];
            model.x = x;
            model.c = x(2)+0.2;
            model.f_x_unconstrained = [x(2);-x(1)];
            model.x0 = [0;pi-2];
        end

        function model = pds_ocp(tc)
            import casadi.*
            
            model = nosnoc.model.Pds();
            x = SX.sym('x',2);
            model.x = [x];
            model.lbx = [-inf;-inf];
            model.ubx = [inf;inf];
            x0 =[1; 5];
            model.x0 = [x0];
            u = SX.sym('u1', 2);;
            model.u = [u];
            model.lbu = [-1;-1];
            model.ubu = [1;1];
            model.u0 = [0;0];
            P = [1/4, 0;
                0, 1/16];
            model.c = [x'*P*x - 1];
            model.f_x_unconstrained = [-0.2*(x(1)+1)^2;-0.4*(x(2)+3)] + u;
        end

        function model = open_gate(tc)
            import casadi.*

            model = nosnoc.model.Pds();
            
            T = 5;
            R = 1;
            R_obj = 1;
            R_obstacle = 5;
            %% Define projected system
            x1 = SX.sym('x1', 2);
            x2 = SX.sym('x2', 2);
            x3 = SX.sym('x3', 2);
            gate = SX.sym('gate', 1);
            x = [x1;x2;x3;gate];
            x_target = [-10;3;-10;0;-7;3;0];
            model.x = x;
            model.lbx = [-inf;-inf;-inf;-inf;-inf;-inf;-inf];
            model.ubx = [inf;5;inf;5;inf;5;inf];
            model.x0 = [-10;2;10;4;7;3.5;0.9];
            u1 = SX.sym('u1', 2);
            u2 = SX.sym('u2', 2);
            model.u = [u1;u2];
            model.lbu = [-10;-10;-10;-10];
            model.ubu = [10;10;10;10];
            model.u0 = [0;0;0;0];
            model.c = [norm_2(x3-x1)-(R+R_obj);
                norm_2(x3-x2)-(R+R_obj);
                norm_2(x2-x1)-(R+R);
                x1(2)-gate-R];
            model.g_path = [x2(2)-gate-R;
                x3(2)-gate-R_obj;
                norm_2(x1-[0;5]) - R-R_obstacle;
                norm_2(x2-[0;5]) - R-R_obstacle;
                norm_2(x3-[0;5]) - R_obj-R_obstacle];
            model.lbg_path = [0;0;0;0;0];
            model.ubg_path = [inf;inf;inf;inf;inf];
            model.f_x_unconstrained = [u1;u2;0;0;0];
        end
    end
end
