classdef TestModelBase < matlab.unittest.TestCase
    methods(Test)
        function test_x_missing(tc)
            opts = nosnoc.Options();
            opts.N_stages = 10;
            model = nosnoc.model.Base();

            tc.verifyError(@()model.verify_and_backfill(opts), 'nosnoc:model:Base:x_missing')
        end

        function test_correct_bound_and_init_population1(tc)
            import casadi.*

            % Test x
            opts1 = nosnoc.Options();
            model1 = nosnoc.model.Base();
            model1.x = SX.sym('x', 3);
            model1.u = SX.sym('u', 3);
            model1.z = SX.sym('z', 3);
            model1.v_global = SX.sym('v', 3);

            model1.verify_and_backfill(opts1);

            % x
            tc.verifySize(model1.lbx, [3,1]);
            tc.verifySize(model1.ubx, [3,1]);
            tc.verifySize(model1.x0, [3,1]);

            % u
            tc.verifySize(model1.lbu, [3,1]);
            tc.verifySize(model1.ubu, [3,1]);
            tc.verifySize(model1.u0, [3,1]);

            % z
            tc.verifySize(model1.lbz, [3,1]);
            tc.verifySize(model1.ubz, [3,1]);
            tc.verifySize(model1.z0, [3,1]);

            % v_global
            tc.verifySize(model1.lbv_global, [3,1]);
            tc.verifySize(model1.ubv_global, [3,1]);
            tc.verifySize(model1.v0_global, [3,1]);
        end
        
        function test_correct_bound_and_init_population2(tc)
            import casadi.*

            % Test x
            opts1 = nosnoc.Options();
            model1 = nosnoc.model.Base();
            model1.x = SX.sym('x', 3);

            model1.verify_and_backfill(opts1);

            % x
            tc.verifySize(model1.lbx, [3,1]);
            tc.verifySize(model1.ubx, [3,1]);
            tc.verifySize(model1.x0, [3,1]);

            % u
            tc.verifyEmpty(model1.lbu);
            tc.verifyEmpty(model1.ubu);
            tc.verifyEmpty(model1.u0);

            % z
            tc.verifyEmpty(model1.lbz);
            tc.verifyEmpty(model1.ubz);
            tc.verifyEmpty(model1.z0);

            % v_global
            tc.verifyEmpty(model1.lbv_global);
            tc.verifyEmpty(model1.ubv_global);
            tc.verifyEmpty(model1.v0_global);
        end

        function test_dims_struct(tc)
            import casadi.*

            % Test x
            opts1 = nosnoc.Options();
            model1 = nosnoc.model.Base();
            model1.x = SX.sym('x', 1);
            model1.u = SX.sym('u', 2);
            model1.z = SX.sym('z', 3);
            model1.v_global = SX.sym('v', 4);

            model1.verify_and_backfill(opts1);

            % x
            tc.assertTrue(isfield(model1.dims, 'n_x'));
            tc.assertTrue(isfield(model1.dims, 'n_u'));
            tc.assertTrue(isfield(model1.dims, 'n_z'));
            tc.assertTrue(isfield(model1.dims, 'n_v_global'));
            tc.assertTrue(isfield(model1.dims, 'n_p_global'));
            tc.assertTrue(isfield(model1.dims, 'n_p_time_var'));
            tc.assertTrue(isfield(model1.dims, 'n_g_path'));
            tc.assertTrue(isfield(model1.dims, 'n_g_comp'));
            tc.assertTrue(isfield(model1.dims, 'n_g_terminal'));

            tc.verifyEqual(model1.dims.n_x, 1);
            tc.verifyEqual(model1.dims.n_u, 2);
            tc.verifyEqual(model1.dims.n_z, 3);
            tc.verifyEqual(model1.dims.n_v_global, 4);
            tc.verifyEqual(model1.dims.n_p_global, 0);
            tc.verifyEqual(model1.dims.n_p_time_var, 0);
            tc.verifyEqual(model1.dims.n_g_path, 0);
            tc.verifyEqual(model1.dims.n_g_comp, 0);
            tc.verifyEqual(model1.dims.n_g_terminal, 0);
            
        end
    end
end
