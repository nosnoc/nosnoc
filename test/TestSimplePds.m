classdef TestSimplePds < matlab.unittest.TestCase
% Watchdog test for PDS->GCS pipeline
    properties (TestParameter)
        rk_representation = {'differential_lift_x','integral'};
        cross_comp_mode = {1, 3, 4, 7}; 
    end
    
    methods (Test)
    end
    methods (Test, ParameterCombination = 'exhaustive')
        function test_fesd_integrator(testCase,rk_representation, cross_comp_mode)
            
            [x_res,t_grid,model,problem_options, solver_options] = test_simple_pds(rk_representation, cross_comp_mode);

            vec = x_res(:,end);
            testCase.verifyLessThan(min(vecnorm(vec-[-1;0], 2)), 1e-5);
        end
    end 
end
