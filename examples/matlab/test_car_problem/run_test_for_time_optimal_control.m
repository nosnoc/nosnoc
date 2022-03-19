function [ ] = run_test_for_time_optimal_control()
%% This function runs the test for different scenarios of time optimal control problems without time-freezing and creates a log file.
%% Open file
fileID = fopen('run_tests.txt','w');
fprintf(fileID,'Running tests on time optimal control problems with FESD (without time-freezing). \n');
fprintf(fileID,'Description of scenarios. \n');
fprintf(fileID,'0: fixed T, just reach final position. \n');
fprintf(fileID,'A: time_optimal_control_problem = 1, equidistant_control_grid = 0, use_speed_of_time_variables = 0, local_speed_of_time_variable = 0. \n');
fprintf(fileID,'B: time_optimal_control_problem = 1, equidistant_control_grid = 1, use_speed_of_time_variables = 0, local_speed_of_time_variable = 0. \n');
fprintf(fileID,'C: time_optimal_control_problem = 1, equidistant_control_grid = 0, use_speed_of_time_variables = 1, local_speed_of_time_variable = 0. \n');
fprintf(fileID,'D: time_optimal_control_problem = 1, equidistant_control_grid = 1, use_speed_of_time_variables = 1, local_speed_of_time_variable = 0. \n');
fprintf(fileID,'E: time_optimal_control_problem = 1, equidistant_control_grid = 0, use_speed_of_time_variables = 1, local_speed_of_time_variable = 1. \n');
fprintf(fileID,'F: time_optimal_control_problem = 1, equidistant_control_grid = 1, use_speed_of_time_variables = 1, local_speed_of_time_variable = 1. \n');
fprintf(fileID,'------------------------\n');
fclose(fileID);
%%
fileID = fopen('run_tests.txt','a+');
fprintf(fileID,'Running tests on time optimal control problems with FESD (without time-freezing). \n');
fprintf(fileID,'------------------------\n');
fclose(fileID);
%% Run Scenarios
% collocation settings


settings.opts_ipopt.ipopt.print_level = 0;

%% Scenario vectors

time_optimal_problem_vec = [0 1 1 1 1 1 1];
equidistant_control_grid_vec = [0 0 1 0 1 0 1];
use_speed_of_time_variables_vec = [0 0 0 1 1 1 1];
local_speed_of_time_variable_vec = [0 0 0 0 0 1 1];


N_exp = 7;
line_str = '-----------------------------------------\n';

% OPEN
fileID = fopen('run_tests.txt','a+');
total_time = 0;
experiments_successful = 0;
tic;
%% Run Experiments
for ii = 1:N_exp


    settings.time_optimal_problem = time_optimal_problem_vec(ii);
    settings.use_speed_of_time_variables =  use_speed_of_time_variables_vec(ii);
    settings.local_speed_of_time_variable = local_speed_of_time_variable_vec(ii);
    settings.equidistant_control_grid = equidistant_control_grid_vec(ii);

    test_message_str = ['Running test ' num2str(ii) ' of ' num2str(N_exp) '.\n'];

    fprintf(test_message_str);
    fprintf(line_str);
    pause(0.5)
    [results,stats] = solve_simple_car_ocp(settings);



    fprintf(fileID,test_message_str);
    fprintf(fileID,line_str);
    if results.status == 1
        fprintf(fileID,'Test sucessful with complementarity residual = %2.2e.\n',stats.complementarity_stats(end));
        fprintf('Current test sucessful!\n');
        experiments_successful = experiments_successful+1;
    else
        fprintf('Current test failed.\n');

    end

    field_names_msg = fieldnames(results.messages);
    N_msg = length(field_names_msg);
    for jj = 1:N_msg;   eval(['fprintf(fileID,results.messages.' field_names_msg{jj} ');']); end

end
total_time = toc;

% Close LOG file
fprintf('---------- SUCESS:  ALL TESTS FINISHED!------------\n' );
fprintf('---------- SUCESS:  ALL TESTS FINISHED!------------\n' );
fprintf('---------- SUCESS:  ALL TESTS FINISHED!------------\n' );
fprintf('Homotopy loop converged in %d of %d experiments. \n',experiments_successful,N_exp);
fprintf('Total experiment time: %2.2f seconds. \n',total_time);

fprintf(fileID,'---------- SUCESS:  ALL TESTS FINISHED!------------ \n');
fprintf(fileID,'Total experiment time: %2.2f seconds. \n',total_time);
fprintf(fileID,'Homotopy loop converged in %d of %d experiments. \n',experiments_successful,N_exp);

fclose(fileID);
end

