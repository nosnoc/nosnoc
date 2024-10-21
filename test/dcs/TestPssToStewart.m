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
        end
    end
end
