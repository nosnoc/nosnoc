classdef TestInelasticMulticontact < matlab.unittest.TestCase
% Testing time freezing for inelastic contacts with friction.
    methods(Test)
        function test_basic_default_opts(tc)
            import casadi.*

            model = tc.simplest_model();
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            heaviside_model = nosnoc.time_freezing.cls_inelastic_multicontact(model, opts);

            tc.verifyEqual(size(heaviside_model.alpha), [2,1]);
            tc.verifyEqual(str(heaviside_model.c), '[q, v]');
            tc.verifyEqual(str(heaviside_model.f_x), '[(v*theta_0), ((-9.81*theta_0)+(100*theta_1)), theta_0]');
            tc.verifyEqual(str(heaviside_model.x), '[q, v, t]');
            tc.verifyEqual(str(heaviside_model.z), '[theta_0, theta_1, beta_bilinear_ode]');
            tc.verifyEqual(str(heaviside_model.g_z), '@1=1, [(theta_0-((alpha_0+alpha_1)-beta_bilinear_ode)), (theta_1-((@1-alpha_0)*(@1-alpha_1))), (beta_bilinear_ode-(alpha_0*alpha_1))]');

            tc.verifyEqual(heaviside_model.x0, [model.x0;0]);
        end

        function test_basic_stabilizing_dynamics(tc)
            import casadi.*

            model = tc.simplest_model();
            opts = nosnoc.Options();
            opts.stabilizing_q_dynamics = true;
            opts.preprocess();
            model.verify_and_backfill(opts);

            heaviside_model = nosnoc.time_freezing.cls_inelastic_multicontact(model, opts);

            tc.verifyEqual(size(heaviside_model.alpha), [2,1]);
            tc.verifyEqual(str(heaviside_model.c), '[q, v]');
            tc.verifyEqual(str(heaviside_model.f_x), '[((v*theta_0)+((-1e-05*q)*theta_1)), ((-9.81*theta_0)+(100*theta_1)), theta_0]');
            tc.verifyEqual(str(heaviside_model.x), '[q, v, t]');
            tc.verifyEqual(str(heaviside_model.z), '[theta_0, theta_1, beta_bilinear_ode]');
            tc.verifyEqual(str(heaviside_model.g_z), '@1=1, [(theta_0-((alpha_0+alpha_1)-beta_bilinear_ode)), (theta_1-((@1-alpha_0)*(@1-alpha_1))), (beta_bilinear_ode-(alpha_0*alpha_1))]');

            tc.verifyEqual(heaviside_model.x0, [model.x0;0]);
        end

        function test_basic_different_M(tc)
            import casadi.*

            model = tc.simplest_model();
            model.M = 3;
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            heaviside_model = nosnoc.time_freezing.cls_inelastic_multicontact(model, opts);

            tc.verifyEqual(size(heaviside_model.alpha), [2,1]);
            tc.verifyEqual(str(heaviside_model.c), '[q, v]');
            tc.verifyEqual(str(heaviside_model.f_x), '[(v*theta_0), ((-3.27*theta_0)+(33.3333*theta_1)), theta_0]');
            tc.verifyEqual(str(heaviside_model.x), '[q, v, t]');
            tc.verifyEqual(str(heaviside_model.z), '[theta_0, theta_1, beta_bilinear_ode]');
            tc.verifyEqual(str(heaviside_model.g_z), '@1=1, [(theta_0-((alpha_0+alpha_1)-beta_bilinear_ode)), (theta_1-((@1-alpha_0)*(@1-alpha_1))), (beta_bilinear_ode-(alpha_0*alpha_1))]');

            tc.verifyEqual(heaviside_model.x0, [model.x0;0]);
        end

        function test_basic_different_a_n(tc)
            import casadi.*

            model = tc.simplest_model();
            opts = nosnoc.Options();
            opts.a_n = 50;
            opts.preprocess();
            model.verify_and_backfill(opts);

            heaviside_model = nosnoc.time_freezing.cls_inelastic_multicontact(model, opts);

            tc.verifyEqual(size(heaviside_model.alpha), [2,1]);
            tc.verifyEqual(str(heaviside_model.c), '[q, v]');
            tc.verifyEqual(str(heaviside_model.f_x), '[(v*theta_0), ((-9.81*theta_0)+(50*theta_1)), theta_0]');
            tc.verifyEqual(str(heaviside_model.x), '[q, v, t]');
            tc.verifyEqual(str(heaviside_model.z), '[theta_0, theta_1, beta_bilinear_ode]');
            tc.verifyEqual(str(heaviside_model.g_z), '@1=1, [(theta_0-((alpha_0+alpha_1)-beta_bilinear_ode)), (theta_1-((@1-alpha_0)*(@1-alpha_1))), (beta_bilinear_ode-(alpha_0*alpha_1))]');

            tc.verifyEqual(heaviside_model.x0, [model.x0;0]);
        end

        function test_basic_quadrature_state(tc)
            import casadi.*

            model = tc.simplest_model();
            model.f_q = model.x(1).^2;
            opts = nosnoc.Options();
            opts.time_freezing_quadrature_state = true;
            opts.preprocess();
            model.verify_and_backfill(opts);

            heaviside_model = nosnoc.time_freezing.cls_inelastic_multicontact(model, opts);

            tc.verifyEqual(size(heaviside_model.alpha), [2,1]);
            tc.verifyEqual(str(heaviside_model.c), '[q, v]');
            tc.verifyEqual(str(heaviside_model.f_x), '[(v*theta_0), ((-9.81*theta_0)+(100*theta_1)), (sq(q)*theta_0), theta_0]');
            tc.verifyEqual(str(heaviside_model.x), '[q, v, L, t]');
            tc.verifyEqual(str(heaviside_model.z), '[theta_0, theta_1, beta_bilinear_ode]');
            tc.verifyEqual(str(heaviside_model.g_z), '@1=1, [(theta_0-((alpha_0+alpha_1)-beta_bilinear_ode)), (theta_1-((@1-alpha_0)*(@1-alpha_1))), (beta_bilinear_ode-(alpha_0*alpha_1))]');

            tc.verifyEqual(heaviside_model.x0, [model.x0;0;0]);
        end

        function test_friction_default_opts(tc)
            import casadi.*

            model = tc.simplest_friction();
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            heaviside_model = nosnoc.time_freezing.cls_inelastic_multicontact(model, opts);

            tc.verifyEqual(size(heaviside_model.alpha), [3,1]);
            tc.verifyEqual(str(heaviside_model.c), '[q_1, v_1, v_0]');
            tc.verifyEqual(str(heaviside_model.f_x), '@1=100, [(v_0*theta_0), (v_1*theta_0), ((-10*theta_1)+(10*theta_2)), (((-9.81*theta_0)+(@1*theta_1))+(@1*theta_2)), theta_0]');
            tc.verifyEqual(str(heaviside_model.x), '[q_0, q_1, v_0, v_1, t]');
            tc.verifyEqual(str(heaviside_model.z), '[theta_0, theta_1, theta_2, beta_bilinear_ode, beta_bilinear_aux]');
            tc.verifyEqual(str(heaviside_model.g_z), '@1=1, [(theta_0-((alpha_0+alpha_1)-beta_bilinear_ode)), (theta_1-(beta_bilinear_aux*alpha_2)), (theta_2-(beta_bilinear_aux*(@1-alpha_2))), (beta_bilinear_ode-(alpha_0*alpha_1)), (beta_bilinear_aux-((@1-alpha_0)*(@1-alpha_1)))]');

            tc.verifyEqual(heaviside_model.x0, [model.x0;0]);
        end

        function test_two_contacts_default_opts(tc)
            import casadi.*

            model = tc.two_contacts();
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            heaviside_model = nosnoc.time_freezing.cls_inelastic_multicontact(model, opts);

            tc.verifyEqual(size(heaviside_model.alpha), [4,1]);
            tc.verifyEqual(str(heaviside_model.c), '[q_0, q_1, v_0, v_1]');
            tc.verifyEqual(str(heaviside_model.f_x), '@1=100, [(v_0*theta_0), (v_1*theta_0), (@1*theta_1), ((-9.81*theta_0)+(@1*theta_2)), theta_0]');
            tc.verifyEqual(str(heaviside_model.x), '[q_0, q_1, v_0, v_1, t]');
            tc.verifyEqual(str(heaviside_model.z), '[theta_0, theta_1, theta_2, beta_bilinear_ode_0, beta_bilinear_ode_1]');
            tc.verifyEqual(str(heaviside_model.g_z), '@1=1, [(theta_0-(((alpha_0+alpha_2)-beta_bilinear_ode_0)*((alpha_1+alpha_3)-beta_bilinear_ode_1))), (theta_1-((@1-alpha_0)*(@1-alpha_2))), (theta_2-((@1-alpha_1)*(@1-alpha_3))), (beta_bilinear_ode_0-(alpha_0*alpha_2)), (beta_bilinear_ode_1-(alpha_1*alpha_3))]');

            tc.verifyEqual(heaviside_model.x0, [model.x0;0]);
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
