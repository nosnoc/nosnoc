
T = {2,2,3,3,2.5};
N_FE = {2,3,2,2,2};
N_stg = {20,20,30,20,25};
max_iter = {1e3,1e3,1e3,1e3,1e3};
cross_comp_mode = {1,1,1,1,3};
q_x_final = {2,2,2,3,2.5};
IrkScheme = {IRKSchemes.RADAU_IIA,IRKSchemes.GAUSS_LEGENDRE,IRKSchemes.RADAU_IIA,IRKSchemes.RADAU_IIA,IRKSchemes.RADAU_IIA};

scenario.Q = diag([50; 10;  10; 10*ones(n_q-3,1); 0.1*ones(n_q,1)]);
scenario.Q_terminal = diag([500; 500; 500; 100*ones(n_q-3,1); 50*ones(n_q,1)]);
scenario.objective = 'walk';
scenario.sigma_0 = 1e2;
scenario.N_homotopy = 8;
scenario.R = diag(0.1*ones(n_u,1));
for ii = 1:length(T)
    scenario.T = T{ii};
    scenario.N_FE = N_FE{ii};
    scenario.N_stg = N_stg{ii};
    scenario.filename = ['quad_half_take' num2str(ii)];
    scenario.max_iter = max_iter{ii};
    scenario.cross_comp_mode = cross_comp_mode{ii};
    scenario.q_x_final = q_x_final{ii};
    scenario.IrkScheme = IrkScheme{ii};
    [model,results] = run_quadroped_half_exp(scenario);
end