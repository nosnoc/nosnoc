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
    end
end
