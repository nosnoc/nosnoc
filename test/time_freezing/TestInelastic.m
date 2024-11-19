classdef TestInelastic < matlab.unittest.TestCase
% Testing time freezing for inelastic contacts with friction.
    methods(Test)
        function test_basic_default_opts(tc)
            import casadi.*

            model = tc.simplest_model();
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            pss_model = nosnoc.time_freezing.cls_inelastic(model, opts);

            tc.verifyEqual(pss_model.S, [1,1;1,-1;-1,1;-1, -1]);
            tc.verifyEqual(str(pss_model.c), '[q, v]');
            tc.verifyEqual(str(pss_model.x), '[q, v, t]');

            F_correct = char(join(["@1=-9.81, @2=1, @3=0, "
                                    "[[v, v, v, @3], "
                                    " [@1, @1, @1, 100], "
                                    " [@2, @2, @2, @3]]"], newline));
            tc.verifyEqual(str(pss_model.F), F_correct);

            tc.verifyEqual(pss_model.x0, [model.x0;0]);
        end

        function test_basic_stabilizing_dynamics(tc)
            import casadi.*

            model = tc.simplest_model();
            opts = nosnoc.Options();
            opts.stabilizing_q_dynamics = true;
            opts.preprocess();
            model.verify_and_backfill(opts);

            pss_model = nosnoc.time_freezing.cls_inelastic(model, opts);

            tc.verifyEqual(pss_model.S, [1,1;1,-1;-1,1;-1, -1]);
            tc.verifyEqual(str(pss_model.c), '[q, v]');
            tc.verifyEqual(str(pss_model.x), '[q, v, t]');

            F_correct = char(join(["@1=-9.81, @2=1, "
                                    "[[v, v, v, (-1e-05*q)], "
                                    " [@1, @1, @1, 100], "
                                    " [@2, @2, @2, 0]]"], newline));
            tc.verifyEqual(str(pss_model.F), F_correct);

            tc.verifyEqual(pss_model.x0, [model.x0;0]);
        end

        function test_basic_different_M(tc)
            import casadi.*

            model = tc.simplest_model();
            model.M = 3;
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            pss_model = nosnoc.time_freezing.cls_inelastic(model, opts);

            tc.verifyEqual(pss_model.S, [1,1;1,-1;-1,1;-1, -1]);
            tc.verifyEqual(str(pss_model.c), '[q, v]');
            tc.verifyEqual(str(pss_model.x), '[q, v, t]');

            F_correct = char(join(["@1=-3.27, @2=1, @3=0, "
                                    "[[v, v, v, @3], "
                                    " [@1, @1, @1, 33.3333], "
                                    " [@2, @2, @2, @3]]"], newline));
            tc.verifyEqual(str(pss_model.F), F_correct);

            tc.verifyEqual(pss_model.x0, [model.x0;0]);
        end

        function test_basic_different_a_n(tc)
            import casadi.*

            model = tc.simplest_model();
            opts = nosnoc.Options();
            opts.a_n = 50;
            opts.preprocess();
            model.verify_and_backfill(opts);

            pss_model = nosnoc.time_freezing.cls_inelastic(model, opts);

            tc.verifyEqual(pss_model.S, [1,1;1,-1;-1,1;-1, -1]);
            tc.verifyEqual(str(pss_model.c), '[q, v]');
            tc.verifyEqual(str(pss_model.x), '[q, v, t]');

            F_correct = char(join(["@1=-9.81, @2=1, @3=0, "
                                    "[[v, v, v, @3], "
                                    " [@1, @1, @1, 50], "
                                    " [@2, @2, @2, @3]]"], newline));
            tc.verifyEqual(str(pss_model.F), F_correct);

            tc.verifyEqual(pss_model.x0, [model.x0;0]);
        end

        function test_basic_quadrature_state(tc)
            import casadi.*

            model = tc.simplest_model();
            model.f_q = model.x(1).^2;
            opts = nosnoc.Options();
            opts.time_freezing_quadrature_state = true;
            opts.preprocess();
            model.verify_and_backfill(opts);

            pss_model = nosnoc.time_freezing.cls_inelastic(model, opts);

            tc.verifyEqual(pss_model.S, [1,1;1,-1;-1,1;-1, -1]);
            tc.verifyEqual(str(pss_model.c), '[q, v]');
            tc.verifyEqual(str(pss_model.x), '[q, v, L, t]');

            F_correct = char(join(["@1=-9.81, @2=sq(q), @3=1, @4=0, "
                                    "[[v, v, v, @4], "
                                    " [@1, @1, @1, 100], "
                                    " [@2, @2, @2, @4], "
                                    " [@3, @3, @3, @4]]"], newline));
            tc.verifyEqual(str(pss_model.F), F_correct);
            tc.verifyEqual(str(pss_model.f_q_T), 'L');

            tc.verifyEqual(pss_model.x0, [model.x0;0;0]);
        end

        function test_friction_default_opts(tc)
            import casadi.*

            model = tc.simplest_friction();
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            pss_model = nosnoc.time_freezing.cls_inelastic(model, opts);

            tc.verifyEqual(pss_model.S, [1     1     1
                                          1     1    -1
                                          1    -1     1
                                          1    -1    -1
                                          -1     1     1
                                          -1     1    -1
                                          -1    -1     1
                                          -1    -1    -1]);
            tc.verifyEqual(str(pss_model.c), '[q_1, v_1, v_0]');
            tc.verifyEqual(str(pss_model.x), '[q_0, q_1, v_0, v_1, t]');

            F_correct = char(join(["@1=0, @2=-9.81, @3=1, @4=100, "
                                    "[[v_0, v_0, v_0, v_0, v_0, v_0, @1, @1], "
                                    " [v_1, v_1, v_1, v_1, v_1, v_1, @1, @1], "
                                    " [@1, @1, @1, @1, @1, @1, -10, 10], "
                                    " [@2, @2, @2, @2, @2, @2, @4, @4], "
                                    " [@3, @3, @3, @3, @3, @3, @1, @1]]"], newline));
            tc.verifyEqual(str(pss_model.F), F_correct);

            tc.verifyEqual(pss_model.x0, [model.x0;0]);
        end

        function test_two_contacts_default_opts(tc)
            import casadi.*

            model = tc.two_contacts();
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            pss_model = nosnoc.time_freezing.cls_inelastic(model, opts);

            tc.verifyEqual(pss_model.S, [1     1     1     1
                                          1     1     1    -1
                                          1     1    -1     1
                                          1     1    -1    -1
                                          1    -1     1     1
                                          1    -1     1    -1
                                          1    -1    -1     1
                                          1    -1    -1    -1
                                          -1     1     1     1
                                          -1     1     1    -1
                                          -1     1    -1     1
                                          -1     1    -1    -1
                                          -1    -1     1     1
                                          -1    -1     1    -1
                                          -1    -1    -1     1
                                          -1    -1    -1    -1]);
            tc.verifyEqual(str(pss_model.c), '[q_0, q_1, v_0, v_1]');
            tc.verifyEqual(str(pss_model.x), '[q_0, q_1, v_0, v_1, t]');

            F_correct = char(join(["@1=0, @2=-9.81, @3=1, @4=100, "
                                    "[[v_0, v_0, v_0, v_0, v_0, @1, v_0, @1, v_0, v_0, @1, @1, v_0, @1, @1, @1, @1], "
                                    " [v_1, v_1, v_1, v_1, v_1, @1, v_1, @1, v_1, v_1, @1, @1, v_1, @1, @1, @1, @1], "
                                    " [@1, @1, @1, @1, @1, @1, @1, @1, @1, @1, @4, @4, @1, @1, @4, @4, @1], "
                                    " [@2, @2, @2, @2, @2, @4, @2, @4, @2, @2, @1, @1, @2, @4, @1, @1, @4], "
                                    " [@3, @3, @3, @3, @3, @1, @3, @1, @3, @3, @1, @1, @3, @1, @1, @1, @1]]"], newline));
            tc.verifyEqual(str(pss_model.F), F_correct);

            tc.verifyEqual(pss_model.x0, [model.x0;0]);
        end
    end

    methods
        function model = simplest_model(tc)
            import casadi.*
            model = nosnoc.model.Cls();

            g = 9.81;
            x0 = [0.8;0];

            q = SX.sym('q',1);
            v = SX.sym('v',1);
            model.M = 1;
            model.x = [q;v];
            model.e = 0;
            model.mu = 0;
            model.x0 = x0;
            model.f_v = -g;
            model.f_c = q;
            model.dims.n_dim_contact = 1;
        end

        function model = simplest_friction(tc)
            import casadi.*
            model = nosnoc.model.Cls();
            
            g = 9.81;

            q = SX.sym('q',2); 
            v = SX.sym('v',2);

            model.x = [q;v]; 
            model.e = 0;
            model.mu = 0.1;
            model.x0 = [0;0.06;3;0]; 
            model.f_v = [0;-g];
            model.f_c = q(2);
            model.J_tangent = [1; 0];
            model.dims.n_dim_contact = 1;
        end

        function model = two_contacts(tc)
            import casadi.*
            model = nosnoc.model.Cls();

            g = 9.81;

            q = SX.sym('q',2); 
            v = SX.sym('v',2);

            model.x = [q;v]; 
            model.e = 0;
            model.mu = 0.0;
            model.x0 = [1;0.06;3;0]; 
            model.f_v = [0;-g];
            model.f_c = q;
            model.dims.n_dim_contact = 1;
        end
    end
end
