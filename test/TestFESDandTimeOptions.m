classdef TestFESDandTimeOptions < matlab.unittest.TestCase
    properties (TestParameter)
        use_fesd = {0,1};
        time_optimal_problem = {0,1};
        equidistant_control_grid = {0,1};
        use_speed_of_time_variables = {0,1};
        local_speed_of_time_variable = {0,1};
    end
    
    methods (Test, ParameterCombination='exhaustive')
        function test_cross_comp_modes(tc,use_fesd,time_optimal_problem,equidistant_control_grid, use_speed_of_time_variables,local_speed_of_time_variable)
            import matlab.unittest.constraints.IssuesNoWarnings;
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            tc.applyFixture(SuppressedWarningsFixture({"nosnoc:homotopy_solver:NLP_infeasible",...
                                                        "nosnoc:NosnocOptions:erroneous_use_speed_of_time_variables",...
                                                        "nosnoc:NosnocOptions:erroneous_local_speed_of_time_variable",...
                                                        "nosnoc:ocp:Solver:terminal_time_for_non_time_optimal_ocp"}));
            [model,problem_options,solver_options,ocp_solver] = test_fesd_and_time_options(use_fesd, time_optimal_problem, equidistant_control_grid, use_speed_of_time_variables, local_speed_of_time_variable);

            % TODO: this should have some better testing.
        end
    end 
end
