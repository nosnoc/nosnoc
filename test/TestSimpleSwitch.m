classdef TestSimpleSwitch < matlab.unittest.TestCase
    properties (TestParameter)
        rk_representation = {'differential','integral'};
        rk_scheme = {RKSchemes.RADAU_IIA,RKSchemes.GAUSS_LEGENDRE};
        dcs_mode = {'Step','Stewart'};
        cross_comp_mode = {1, 3, 4, 7}; 
    end
    
    methods (Test)
    end
    methods (Test, ParameterCombination = 'exhaustive')
        function test_fesd_integrator(testCase,rk_representation, rk_scheme, dcs_mode, cross_comp_mode)
            import matlab.unittest.constraints.IssuesNoWarnings;
            issuesNoWarningsConstraint = IssuesNoWarnings('WhenNargoutIs', 5);
            testCase.verifyThat(@() test_simple_switch(rk_representation, rk_scheme, dcs_mode, cross_comp_mode), issuesNoWarningsConstraint);

            [results,stats,model,problem_options, solver_options] = issuesNoWarningsConstraint.FunctionOutputs{:};

            vec = [results.x; results.t_grid'];
            testCase.assertLessThan(min(vecnorm(vec-[0;0.5], 2)), 1e-5);
        end
    end 
end
