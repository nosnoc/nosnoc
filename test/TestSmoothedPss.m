classdef TestSmoothedPss < matlab.unittest.TestCase
    properties (TestParameter)
        odesolver = {'ode45', 'ode23', 'ode113', 'ode78', 'ode89', 'ode15s', 'ode23s', 'ode23t', 'ode23tb', 'cvodesnonstiff', 'cvodesstiff', 'idas'}
        % The rest are not tested because they are _extremely_ slow, on the order of minutes.
        % 

        model = {'switch', 'sliding'}
    end

    methods (Test, ParameterCombination = 'exhaustive')
        function test_smoothed_pss_integrator(tc, odesolver, model)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            tc.applyFixture(SuppressedWarningsFixture({'nosnoc:integrator:Options:likely_bad_ode_solver'}));
            tc.assumeFalse(ismember(odesolver, {'ode45', 'ode23', 'ode113', 'ode78', 'ode89', 'cvodesstiff', 'idas'}) && strcmp(model, 'sliding')); % Filter bad performing solvers
            tc.assumeFalse((string(version('-release')) < "2024a") && ismember(odesolver, {'cvodesnonstiff', 'cvodesstiff', 'idas'}))
            
            
            fprintf(['using ' odesolver ' to solve the ' model ' model\n'])
            [x_res,x_star] = test_smoothed_pss(odesolver,model);
            
            tc.verifyLessThan(norm(x_res(:,end) - x_star), 1e-3);
        end
    end 
end
