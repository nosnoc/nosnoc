% n_q = 11;
% n_u = 8;
% scenario.T = 2;
% scenario.N_FE = 3;
% scenario.N_stg = 20;
% scenario.filename = 'quad_full1';
% scenario.objective = 'walk';
% scenario.max_iter = 600;
% scenario.q_x_final = 1.5;
% % scenario.q_z_final = 0.7;
% % scenario.impose_terminal_constraint = 1;
% scenario.Q = diag([50; 10;  10; 10*ones(n_q-3,1); 0.1*ones(n_q,1)]);
% scenario.Q_terminal = diag([500; 500; 500; 100*ones(n_q-3,1); 50*ones(n_q,1)]);
% scenario.R = diag(0.1*ones(n_u,1)); 
% error('TODO: create a run_quadroped_full_exp funtion out of quadroped_full_ocp')
% [model,results] = run_quadroped_full_exp(scenario);

%% take 1
scenario.T = 2;
scenario.N_FE = 2;
scenario.N_stg = 20;
scenario.filename = 'quad_half1';
scenario.objective = 'walk';
scenario.max_iter = 1e3;
scenario.sigma_0 = 1e2;
scenario.N_homotopy = 8;
scenario.cross_comp_mode = 1;
scenario.q_x_final = 1.5;
scenario.Q = diag([50; 10;  10; 10*ones(n_q-3,1); 0.1*ones(n_q,1)]);
scenario.Q_terminal = diag([500; 500; 500; 100*ones(n_q-3,1); 50*ones(n_q,1)]);
scenario.R = diag(0.1*ones(n_u,1)); 
scenario.IrkScheme = IRKSchemes.RADAU_IIA;
[model,results] = run_quadroped_half_exp(scenario);

%% take 2
scenario.T = 2;
scenario.N_FE = 2;
scenario.N_stg = 30;
scenario.filename = 'quad_half2';
scenario.objective = 'walk';
scenario.max_iter = 1e3;
scenario.sigma_0 = 1e1;
scenario.N_homotopy = 7;
scenario.cross_comp_mode = 1;
scenario.q_x_final = 1.5;
scenario.Q = diag([50; 10;  10; 10*ones(n_q-3,1); 0.1*ones(n_q,1)]);
scenario.Q_terminal = diag([500; 500; 500; 100*ones(n_q-3,1); 50*ones(n_q,1)]);
scenario.R = diag(0.1*ones(n_u,1)); 
scenario.IrkScheme = IRKSchemes.RADAU_IIA;
[model,results] = run_quadroped_half_exp(scenario);

%% take 3
scenario.T = 2.5;
scenario.N_FE = 2;
scenario.N_stg = 25;
scenario.filename = 'quad_half3';
scenario.objective = 'walk';
scenario.max_iter = 1e3;
scenario.sigma_0 = 1e2;
scenario.N_homotopy = 8;
scenario.cross_comp_mode = 1;
scenario.q_x_final = 2;
scenario.Q = diag([50; 10;  10; 10*ones(n_q-3,1); 0.1*ones(n_q,1)]);
scenario.Q_terminal = diag([500; 500; 500; 100*ones(n_q-3,1); 50*ones(n_q,1)]);
scenario.R = diag(0.1*ones(n_u,1)); 
scenario.IrkScheme = IRKSchemes.RADAU_IIA;
[model,results] = run_quadroped_half_exp(scenario);
%% take 4
scenario.T = 2.5;
scenario.N_FE = 2;
scenario.N_stg = 25;
scenario.filename = 'quad_half4';
scenario.objective = 'walk';
scenario.max_iter = 500;
scenario.sigma_0 = 1e2;
scenario.N_homotopy = 8;
scenario.cross_comp_mode = 1;
scenario.q_x_final = 2;
scenario.Q = diag([50; 10;  10; 10*ones(n_q-3,1); 0.1*ones(n_q,1)]);
scenario.Q_terminal = diag([500; 500; 500; 100*ones(n_q-3,1); 50*ones(n_q,1)]);
scenario.R = diag(0.1*ones(n_u,1)); 
scenario.IrkScheme = IRKSchemes.RADAU_IIA;
[model,results] = run_quadroped_half_exp(scenario);

%% take 5
scenario.T = 2.5;
scenario.N_FE = 2;
scenario.N_stg = 25;
scenario.filename = 'quad_half5';
scenario.objective = 'walk';
scenario.max_iter = 1e3;
scenario.sigma_0 = 1e2;
scenario.N_homotopy = 8;
scenario.cross_comp_mode = 3;
scenario.q_x_final = 2;
scenario.Q = diag([50; 10;  10; 10*ones(n_q-3,1); 0.1*ones(n_q,1)]);
scenario.Q_terminal = diag([500; 500; 500; 100*ones(n_q-3,1); 50*ones(n_q,1)]);
scenario.R = diag(0.1*ones(n_u,1)); 
scenario.IrkScheme = IRKSchemes.RADAU_IIA;
[model,results] = run_quadroped_half_exp(scenario);

%% take 6
scenario.T = 2.5;
scenario.N_FE = 2;
scenario.N_stg = 25;
scenario.filename = 'quad_half6';
scenario.objective = 'walk';
scenario.max_iter = 500;
scenario.sigma_0 = 1e2;
scenario.N_homotopy = 8;
scenario.cross_comp_mode = 1;
scenario.q_x_final = 2;
scenario.Q = diag([50; 10;  10; 10*ones(n_q-3,1); 0.1*ones(n_q,1)]);
scenario.Q_terminal = diag([500; 500; 500; 100*ones(n_q-3,1); 50*ones(n_q,1)]);
scenario.R = diag(0.1*ones(n_u,1)); 
scenario.IrkScheme = IRKSchemes.GAUSS_LEGENDRE;
[model,results] = run_quadroped_half_exp(scenario);

%% take 7
scenario.T = 2.5;
scenario.N_FE = 3;
scenario.N_stg = 25;
scenario.filename = 'quad_half7';
scenario.objective = 'walk';
scenario.max_iter = 500;
scenario.sigma_0 = 1e2;
scenario.N_homotopy = 8;
scenario.cross_comp_mode = 1;
scenario.q_x_final = 2;
scenario.Q = diag([50; 10;  10; 10*ones(n_q-3,1); 0.1*ones(n_q,1)]);
scenario.Q_terminal = diag([500; 500; 500; 100*ones(n_q-3,1); 50*ones(n_q,1)]);
scenario.R = diag(0.1*ones(n_u,1)); 
scenario.IrkScheme = IRKSchemes.RADAU_IIA;
[model,results] = run_quadroped_half_exp(scenario);

%% take 8
scenario.T = 3;
scenario.N_FE = 2;
scenario.N_stg = 30;
scenario.filename = 'quad_half8';
scenario.objective = 'walk';
scenario.max_iter = 1000;
scenario.sigma_0 = 1e2;
scenario.N_homotopy = 8;
scenario.cross_comp_mode = 1;
scenario.q_x_final = 3;
scenario.Q = diag([50; 10;  10; 10*ones(n_q-3,1); 0.1*ones(n_q,1)]);
scenario.Q_terminal = diag([500; 500; 500; 100*ones(n_q-3,1); 50*ones(n_q,1)]);
scenario.R = diag(0.1*ones(n_u,1)); 
scenario.IrkScheme = IRKSchemes.RADAU_IIA;
[model,results] = run_quadroped_half_exp(scenario);

%% nicer implement
quadroped_looper