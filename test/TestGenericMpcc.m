classdef TestGenericMpcc < matlab.unittest.TestCase
    properties (TestParameter)
        mpcc_method = cellstr(enumeration('nosnoc.solver.MpccMethod'));
        elasticity_mode = cellstr(enumeration('ElasticityMode'))
    end
    
    methods (Test)
    end
    methods (Test, ParameterCombination = 'exhaustive')
        function test_fesd_integrator(testCase, mpcc_method, elasticity_mode)
            import matlab.unittest.constraints.IssuesNoWarnings;
            issuesNoWarningsConstraint = IssuesNoWarnings('WhenNargoutIs', 2);
            testCase.verifyThat(@() test_generic_mpcc(mpcc_method, elasticity_mode), issuesNoWarningsConstraint);

            [results,stats] = issuesNoWarningsConstraint.FunctionOutputs{:};

            results.x
            testCase.assertLessThan(vecnorm(results.x-[1;0]), 1e-3);
        end
    end 
end
