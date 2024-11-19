classdef TestElastic < matlab.unittest.TestCase
    methods(Test)
        function test_basic_default_opts(tc)
            import casadi.*

            model = tc.simplest_model();
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            pss_model = nosnoc.time_freezing.cls_elastic(model, opts);

            tc.verifyEqual(pss_model.S, [1;-1]);
            tc.verifyEqual(str(pss_model.c), 'q');
            tc.verifyEqual(str(pss_model.x), '[q, v, t]');

            F_correct = char(join([""
                                    "[[v, v], "
                                    " [-9.81, (-10*q)], "
                                    " [1, 0]]"], newline));
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

            pss_model = nosnoc.time_freezing.cls_elastic(model, opts);

            tc.verifyEqual(pss_model.S, [1;-1]);
            tc.verifyEqual(str(pss_model.c), 'q');
            tc.verifyEqual(str(pss_model.x), '[q, v, t]');

            F_correct = char(join([""
                                    "[[v, (0.333333*v)], "
                                    " [-3.27, (-3.33333*q)], "
                                    " [1, 0]]"], newline));
            tc.verifyEqual(str(pss_model.F), F_correct);

            tc.verifyEqual(pss_model.x0, [model.x0;0]);
        end

        function test_basic_different_k_aux(tc)
            import casadi.*

            model = tc.simplest_model();
            opts = nosnoc.Options();
            opts.k_aux = 5;
            opts.preprocess();
            model.verify_and_backfill(opts);

            pss_model = nosnoc.time_freezing.cls_elastic(model, opts);

            tc.verifyEqual(pss_model.S, [1;-1]);
            tc.verifyEqual(str(pss_model.c), 'q');
            tc.verifyEqual(str(pss_model.x), '[q, v, t]');

            F_correct = char(join([""
                                    "[[v, v], "
                                    " [-9.81, (-5*q)], "
                                    " [1, 0]]"], newline));
            tc.verifyEqual(str(pss_model.F), F_correct);

            tc.verifyEqual(pss_model.x0, [model.x0;0]);
        end

        function test_basic_different_e(tc)
            import casadi.*

            model = tc.simplest_model();
            model.e = 0.5;
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            pss_model = nosnoc.time_freezing.cls_elastic(model, opts);

            tc.verifyEqual(pss_model.S, [1;-1]);
            tc.verifyEqual(str(pss_model.c), 'q');
            tc.verifyEqual(str(pss_model.x), '[q, v, t]');

            F_correct = char(join([""
                                    "[[v, v], "
                                    " [-9.81, ((-10*q)+(-1.36265*v))], "
                                    " [1, 0]]"], newline));
            tc.verifyEqual(str(pss_model.F), F_correct);

            tc.verifyEqual(pss_model.x0, [model.x0;0]);
        end

        function test_basic_quadrature_state(tc)
            import casadi.*

            model = tc.simplest_model();
            model.f_q = model.q.^2;
            opts = nosnoc.Options();
            opts.time_freezing_quadrature_state = true;
            opts.preprocess();
            model.verify_and_backfill(opts);

            pss_model = nosnoc.time_freezing.cls_elastic(model, opts);

            tc.verifyEqual(pss_model.S, [1;-1]);
            tc.verifyEqual(str(pss_model.c), 'q');
            tc.verifyEqual(str(pss_model.x), '[q, v, L, t]');

            F_correct = char(join(["@1=0, "
                                    "[[v, v], "
                                    " [-9.81, (-10*q)], "
                                    " [@1, @1], "
                                    " [1, @1]]"], newline));
            tc.verifyEqual(str(pss_model.F), F_correct);

            tc.verifyEqual(pss_model.x0, [model.x0;0;0]);
        end

        function test_too_many_contacts(tc)
            import casadi.*

            model = tc.simplest_model();
            model.f_c = [model.f_c; model.x(1)+10];
            opts = nosnoc.Options();
            opts.time_freezing_quadrature_state = true;
            opts.preprocess();
            model.verify_and_backfill(opts);

            tc.verifyError(@() nosnoc.time_freezing.cls_elastic(model, opts), "nosnoc:time_freezing:cls_elastic:too_many_contacts");
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
            model.e = 1;
            model.mu = 0;
            model.x0 = x0;
            model.f_v = -g;
            model.f_c = q;
        end
    end
end
