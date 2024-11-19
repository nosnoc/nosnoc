classdef TestIntegrator < matlab.unittest.TestCase
    properties (TestParameter)
        use_fesd = {0,1};
        rk_representation = {'differential','integral'};
        rk_scheme = {RKSchemes.RADAU_IIA,RKSchemes.GAUSS_LEGENDRE};
        dcs_mode = {'Heaviside','Stewart'};
    end

    methods (Test)
    end
    methods (Test, ParameterCombination = 'exhaustive')
        function test_fesd_integrator(testCase,use_fesd, rk_representation, rk_scheme, dcs_mode)
            test_integrator(use_fesd, rk_representation, rk_scheme, dcs_mode);
        end
    end 
end


