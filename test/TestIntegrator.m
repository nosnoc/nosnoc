classdef TestIntegrator < matlab.unittest.TestCase
    %TESTSETTINGS Summary of this class goes here
    %   Detailed explanation goes here
    properties (TestParameter)
        use_fesd = {0,1};
        rk_representation = {'differential','integral'};
        rk_scheme = {RKSchemes.RADAU_IIA,RKSchemes.GAUSS_LEGENDRE};
        dcs_mode = {'Step','Stewart'};
    end

    methods (Test)
    end
    methods (Test, ParameterCombination = 'exhaustive')
        function test_fesd_integrator(testCase,use_fesd, rk_representation, rk_scheme, dcs_mode)
            import matlab.unittest.constraints.IssuesNoWarnings;
            issuesNoWarningsConstraint = IssuesNoWarnings('WhenNargoutIs', 4);
            testCase.verifyThat(@() test_integrator(use_fesd, rk_representation, rk_scheme, dcs_mode), issuesNoWarningsConstraint);
        end
    end 
end


