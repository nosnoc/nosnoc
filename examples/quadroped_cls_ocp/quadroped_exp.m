n_q = 11;
n_u = 8;
scenario.T = 2;
scenario.N_FE = 3;
scenario.N_stg = 20;
scenario.filename = 'quad_full1';
scenario.objective = 'walk';
scenario.max_iter = 600;
scenario.q_x_final = 1.5;
% scenario.q_z_final = 0.7;
% scenario.impose_terminal_constraint = 1;
scenario.Q = diag([50; 10;  10; 10*ones(n_q-3,1); 0.1*ones(n_q,1)]);
scenario.Q_terminal = diag([500; 500; 500; 100*ones(n_q-3,1); 50*ones(n_q,1)]);
scenario.R = diag(0.1*ones(n_u,1)); 
error('TODO: create a run_quadroped_full_exp funtion out of quadroped_full_ocp')
[model,results] = run_quadroped_full_exp(scenario);

%%
n_q = 7;
n_u = 4;
scenario.a_n = 200;
scenario.T = 2;
scenario.N_FE = 3;
scenario.N_stg = 20;
scenario.filename = 'quad_half1';
scenario.objective = 'walk';
scenario.max_iter = 600;
scenario.q_x_final = 1.5;
% scenario.q_z_final = 0.7;
% scenario.impose_terminal_constraint = 1;
scenario.Q = diag([50; 10;  10; 10*ones(n_q-3,1); 0.1*ones(n_q,1)]);
scenario.Q_terminal = diag([500; 500; 500; 100*ones(n_q-3,1); 50*ones(n_q,1)]);
scenario.R = diag(0.1*ones(n_u,1)); 
[model,results] = run_quadroped_half_exp(scenario);
