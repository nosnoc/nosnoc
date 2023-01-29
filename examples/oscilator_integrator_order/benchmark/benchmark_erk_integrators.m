clc
%% Run experiments to test FESD integrator order for different IRK schemes
N_start = 11;
N_end = 1001;
N_samples = 5;
% N_end = 301;
% N_samples = 3;
%% Explicit RK 
% with fesd
%  
legend_str = {'Explicit Euler', 'Heun 2','Heun 3','Runge-Kutta 4','no','Nystrom 5'};
settings.use_fesd = 1;
settings.irk_scheme = 'Explicit-RK';
settings.irk_representation = 'differential';
settings.save_results = 1;
settings.n_s_vec = [1 2 3 4];
% settings.n_s_vec = [3];
settings.scenario_name = 'erk_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

% without FESD
settings.use_fesd = 0;
settings.scenario_name = 'erk_std_differential';
[results] = integrator_order_experiment(settings,legend_str);
%  

% end