classdef TestGenericMpcc < matlab.unittest.TestCase
    properties (TestParameter)
        mpcc_method = cellstr(enumeration('nosnoc.solver.MpccMethod'));
        homotopy_steering_strategy = cellstr(enumeration('HomotopySteeringStrategy'))
    end
    
    methods (Test)
    end
    methods (Test, ParameterCombination = 'exhaustive')
        function test_fesd_integrator(testCase, mpcc_method, homotopy_steering_strategy)
            import matlab.unittest.constraints.IssuesNoWarnings;
            issuesNoWarningsConstraint = IssuesNoWarnings('WhenNargoutIs', 2);
            testCase.verifyThat(@() test_generic_mpcc(mpcc_method, homotopy_steering_strategy), issuesNoWarningsConstraint);

            [results,stats] = issuesNoWarningsConstraint.FunctionOutputs{:};

            testCase.assertLessThan(vecnorm(results.x-[1;0]), 1e-3);
        end
    end 
end
