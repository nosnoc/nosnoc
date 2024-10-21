classdef TestPssToStewart < matlab.unittest.TestCase
    methods(Test)
        function test_dims_1(testCase)
            import casadi.*
            model = nosnoc.model.Pss();
            model.x = SX.sym('x', 2);
            model.c = x;
            model.S = [1,1;
                -1,1;
                1,-1;
                -1,-1];
        end
    end
end
