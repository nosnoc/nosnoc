classdef TestSettings < matlab.unittest.TestCase
    %TESTSETTINGS Summary of this class goes here
    %   Detailed explanation goes here
    properties (TestParameter)
        cross_comp = {1,2,3};
        mpcc_mode = {1,2,3,4};
    end
    
    methods (Test)
    end
    methods (Test, ParameterCombination = 'pairwise')
        function test_cross_comp_modes(testCase,cross_comp,mpcc_mode)
            import matlab.unittest.constraints.IssuesNoWarnings;
            issuesNoWarningsConstraint = IssuesNoWarnings('WhenNargoutIs', 4);
            testCase.verifyThat(@() test_cross_comp_modes(cross_comp,mpcc_mode), issuesNoWarningsConstraint);

        end
    end 
end


