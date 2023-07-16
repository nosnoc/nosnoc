classdef TestFESDandTimeOptions < matlab.unittest.TestCase
    %TESTSETTINGS Summary of this class goes here
    %   Detailed explanation goes here
    properties (TestParameter)
        use_fesd = {0,1};
        time_optimal_problem = {0,1};
        equidistant_control_grid = {0,1};
        use_speed_of_time_variables = {0,1};
        local_speed_of_time_variable = {0,1};
    end
    
    methods (Test, ParameterCombination='exhaustive')
        function test_cross_comp_modes(testCase,use_fesd,time_optimal_problem,equidistant_control_grid, use_speed_of_time_variables,local_speed_of_time_variable)
            import matlab.unittest.constraints.IssuesNoWarnings;
            warning('off', 'nosnoc:homotopy_solver:NLP_infeasible')
            % Expected warning due to nonsensical settings being handled
            warning('off', 'nosnoc:NosnocOptions:erroneous_use_speed_of_time_variables')
            warning('off', 'nosnoc:NosnocOptions:erroneous_local_speed_of_time_variable')
            issuesNoWarningsConstraint = IssuesNoWarnings('WhenNargoutIs', 5);
            testCase.verifyThat(@() test_fesd_and_time_options(use_fesd, time_optimal_problem, equidistant_control_grid, use_speed_of_time_variables, local_speed_of_time_variable), issuesNoWarningsConstraint);
        end
    end 
end
