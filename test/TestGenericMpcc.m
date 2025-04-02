classdef TestGenericMpcc < matlab.unittest.TestCase
    properties (TestParameter)
        mpcc_method = cellstr(enumeration('nosnoc.reg_homotopy.MpccMethod'));
        homotopy_steering_strategy = cellstr(enumeration('HomotopySteeringStrategy'))
    end
    
    methods (Test)
    end
    methods (Test, ParameterCombination = 'exhaustive')
        function test_fesd_integrator(testCase, mpcc_method, homotopy_steering_strategy)
            [results,stats] = test_generic_mpcc(mpcc_method, homotopy_steering_strategy);

            testCase.assertLessThan(vecnorm(results.x-[1;0]), 1e-3);
        end
    end 
end
