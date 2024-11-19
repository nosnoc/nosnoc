classdef TestSimpleCls < matlab.unittest.TestCase
% Watchdog test for CLS->CLS pipeline
    properties (TestParameter)
        % TODO: differential_lift_x seems very unstable at impacts.
        rk_representation = {'differential','integral'};
        cross_comp_mode = {1, 3, 4, 7}; 
    end
    
    methods (Test)
    end
    methods (Test, ParameterCombination = 'exhaustive')
        function test_cls_integrator(testCase,rk_representation, cross_comp_mode)
            
            [x_res,t_grid,model,problem_options, solver_options] = test_simple_cls(rk_representation, cross_comp_mode);

            vec = x_res(:,end);
            testCase.verifyLessThan(min(vecnorm(vec-[0.5386;2.2645], 2)), 1e-3);
        end
    end 
end
