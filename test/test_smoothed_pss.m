function [x_res,x_star] = test_smoothed_pss(odesolver, model)
    import casadi.*
    problem_options = nosnoc.Options();
    integrator_options = nosnoc.integrator.Options();
    solver_options = integrator_options.fesd_solver_opts;
    
    switch(model)
      case 'switch'
        [x_star, T_sim, N_sim, model] = simple_switch();
      case 'sliding'
        [x_star, T_sim, N_sim, model] = sliding_mode();
    end

    problem_options.N_sim = N_sim;
    problem_options.T_sim = T_sim;

    integrator_options.integrator_plugin = IntegratorType.SMOOTHED_PSS;
    integrator_options.matlab_ode_solver = odesolver;
    
    integrator = nosnoc.Integrator(model, problem_options, integrator_options);
    [t_grid, x_res] = integrator.simulate();
end


function [x_star, T_sim, N_sim, model] = simple_switch()
    import casadi.*
    model = nosnoc.model.Pss();

    N_sim = 16;
    T_sim = 1.5;

    model.x0 = -1;
    x = SX.sym('x',1);
    model.x = x;
    model.c = x;
    model.S = [-1; 1];
    f_1 = [2]; f_2 = [0.2];
    model.F = [f_1 f_2];
    x_star = [0.2];
end

function [x_star, T_sim, N_sim, model] = sliding_mode()
    import casadi.*
    model = nosnoc.model.Pss();
    
    model.x0 = -sqrt(2);
    x = SX.sym('x',1);
    model.x = x;
    model.c = x;
    model.S = [-1; 1];
    f_1 = [1]; f_2 = [-1];
    model.F = [f_1 f_2];

    N_sim = 7;
    T_sim = 2;

    x_star = 0;
end
