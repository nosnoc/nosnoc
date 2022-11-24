classdef TestIntegrator < matlab.unittest.TestCase
    %TESTSETTINGS Summary of this class goes here
    %   Detailed explanation goes here
    properties (TestParameter)
        use_fesd = {0,1};
        irk_representation = {'differential','integral'};
        irk_scheme = {'Radau-IIA','Gauss-Legendre'};
        pss_mode = {'Step','Stewart'};
    end
    
    methods (Test)
    end
    methods (Test, ParameterCombination = 'exhaustive')
        function test_fesd_integrator(testCase,use_fesd, irk_representation, irk_scheme, pss_mode)
            import matlab.unittest.constraints.IssuesNoWarnings;
            issuesNoWarningsConstraint = IssuesNoWarnings('WhenNargoutIs', 4);
            testCase.verifyThat(@() test_integrator(use_fesd, irk_representation, irk_scheme, pss_mode), issuesNoWarningsConstraint);
        end
    end 
end


