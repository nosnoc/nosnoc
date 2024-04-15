classdef TestMpccCrossComp < matlab.unittest.TestCase
    properties (TestParameter)
        cross_comp = {1,3,4,7};
        mpcc_mode = {'direct','Scholtes_eq','Scholtes_ineq','ell_1_penalty',...
             'elastic_ineq','elastic_eq','elastic_two_sided',...
             };
        % NOTE: 'elastic_two_sided' currently fails.
    end
    
    methods (Test, ParameterCombination = 'pairwise')
        function test_cross_comp_modes(testCase,cross_comp,mpcc_mode)
            import matlab.unittest.constraints.IssuesNoWarnings;
            issuesNoWarningsConstraint = IssuesNoWarnings('WhenNargoutIs', 4);
            %testCase.verifyThat(@() test_simple_car_model(cross_comp,mpcc_mode), issuesNoWarningsConstraint);
        end
    end 
end


