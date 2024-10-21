classdef TestSlidingMode < matlab.unittest.TestCase
    properties (TestParameter)
        rk_representation = {'differential','integral'};
        rk_scheme = {RKSchemes.RADAU_IIA,RKSchemes.GAUSS_LEGENDRE};
        dcs_mode = {'Heaviside','Stewart'};
        cross_comp_mode = {1, 3, 4, 7}; 
    end
    
    methods (Test)
    end
    methods (Test, ParameterCombination = 'exhaustive')
        function test_sliding(testCase,rk_representation, rk_scheme, dcs_mode, cross_comp_mode)
            import matlab.unittest.constraints.IssuesNoWarnings;
            issuesNoWarningsConstraint = IssuesNoWarnings('WhenNargoutIs', 6);
            testCase.verifyThat(@() test_sliding_mode(rk_representation, rk_scheme, dcs_mode, cross_comp_mode), issuesNoWarningsConstraint);

            [x_res,t_grid,integrator,model,problem_options,solver_options] = issuesNoWarningsConstraint.FunctionOutputs{:};

            vec = [x_res; t_grid];
            testCase.verifyLessThan(min(vecnorm(vec-[0;sqrt(2)], 2)), 1e-5);
        end
    end 
end
