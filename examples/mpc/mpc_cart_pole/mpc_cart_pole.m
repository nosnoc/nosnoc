clear; clc; close all;
import casadi.*
%% plot and video settings
latexify_plot()
nice_plot_colors;

% visualization
create_video = 0;
plot_intermediate_solutions = 1;
save_qpecs = 0;

linewidth = 1.5;

%% Model settings
% Discretization options
model_plant_missmatch = 1;
model_plant_same_integrator = 0;
n_s = 1;
n_s_sim = 2;
N_sim = 4;

N_FE = 2;
use_fesd = 1;
gamma_h = 1.0;

%% cost weights
Q = diag([20; 10; 1; 1]);
Q = diag([10; 10; 1; 1]);
Q_terminal = Q*100;
R = 0.1;

Q = diag([1; 10; 1; 1]);
Q_terminal = diag([1000; 1000; 1; 1]);
R = 1e-1;
save_data = 1;
%
%% Basics
filename_pdf = 'cartpole_mpc';
scenarios= {'Ideal-MPC','HyRTI','HyAS-RTI-QPCC','HyRTI-IPOPT','HyAS-RTI-QPCC-IPOPT','Smoothing','HyRTI-fast','HyAS-RTI-QP'};
%scenarios= {'Ideal-MPC','HyRTI','Smoothing','HyAS-RTI-QPCC'};

use_ccopt = false; % Recommended: use our fastest QPCC/MPCC solver; requires installation (see nosnoc README).

% for kk = 1:length(scenarios)
for kk = [2]
    scenario =  scenarios{kk};
    % initalize opts
    mpc_options = nosnoc.mpc.Options(); % Initialize all options related to the MPC implementation.
    homotopy_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosonc problems.
    ccopt_options = nosnoc.ccopt.Options(); % ccopt options 
    % MPEC solver options
    tol = 1e-8;
    N_homotopy = 10;
    sigma_0 = 1;
    fast_sigma_0 = 1e-3; % sigma for mpec warmstart in mpc


    % mpecopt_options = mpecopt.Options();
    % mpecopt_options.rho_TR_phase_i_init = 1e-4;
    % mpecopt_options.rho_TR_phase_ii_init = 1e-6;
    % % mpecopt_options.stop_if_S_stationary = true;
    % mpecopt_options.relax_and_project_kappa = 0.1;
    % mpecopt_options.settings_lpec.lpec_solver = "Gurobi";


    switch scenario
        case 'Ideal-MPC'
            data_name = scenario;
            use_rtmpc = 0;
            do_shift_initialization = 0;
            solve_advanced_problem = 0;
            advanced_n_qpecs = 1;
            use_feedback_qp = false;
            use_probing_qp = false;
            warmstart_qpec = false;
            tol = 1e-9;
        case 'HyRTI'
            data_name = scenario;
            use_rtmpc = 1;
            use_ccopt = 0;
            do_shift_initialization = 0;
            solve_advanced_problem = 0;
            advanced_n_qpecs = 1;
            use_feedback_qp = false;
            use_probing_qp = false;
            warmstart_qpec = true;
            % N_homotopy = 5;
            % sigma_0 = 1;
            tol = 1e-8;
        case 'HyRTI-IPOPT'
            data_name = scenario;
            use_rtmpc = 1;
            use_ccopt = false;
            do_shift_initialization = 0;
            solve_advanced_problem = 0;
            advanced_n_qpecs = 1;
            use_feedback_qp = false;
            use_probing_qp = false;
            warmstart_qpec = true;
            % N_homotopy = 5;
            % sigma_0 = 1;
            tol = 1e-8;
        case 'HyRTI-fast'
            data_name = scenario;
            use_rtmpc = 1;
            do_shift_initialization = 0;
            solve_advanced_problem = 0;
            advanced_n_qpecs = 1;
            use_feedback_qp = false;
            use_probing_qp = false;
            warmstart_qpec = 0;
            N_homotopy = 5;
            sigma_0 = 1e-1;
            tol = 1e-4;
        case 'HyAS-RTI-QPCC'
            data_name = scenario;
            use_rtmpc = 1;
            use_ccopt = true;
            do_shift_initialization = 0;
            solve_advanced_problem = 1;
            advanced_n_qpecs = 1;
            use_feedback_qp = false;
            use_probing_qp = false;
            warmstart_qpec = true;
        case 'HyAS-RTI-QPCC-IPOPT'
            data_name = scenario;
            use_rtmpc = 1;
            use_ccopt = false;
            do_shift_initialization = 0;
            solve_advanced_problem = 1;
            advanced_n_qpecs = 1;
            use_feedback_qp = false;
            use_probing_qp = false;
            warmstart_qpec = true;
        case 'HyAS-RTI-QP'
            data_name = scenario;
            use_rtmpc = 1;
            do_shift_initialization = 0;
            solve_advanced_problem = 1;
            advanced_n_qpecs = 1;
            use_feedback_qp = true;
            use_probing_qp = false;
            warmstart_qpec = true;
        case 'Smoothing'
            data_name = scenario;
            use_rtmpc = 0;
            do_shift_initialization = 0;
            solve_advanced_problem = 0;
            advanced_n_qpecs = 1;
            use_feedback_qp = false;
            use_probing_qp = false;
            warmstart_qpec = true;
            N_homotopy = 1;
            sigma_0 = 1e-1;
            tol = 1e-6;
        otherwise
            error('Unknown scenario: %s', scenario);
    end

    % Homotopy options
    homotopy_options.complementarity_tol = tol; % Value to drive the complementarity residual to.
    homotopy_options.N_homotopy = N_homotopy; % Maximum number of homotopy iterations.
    homotopy_options.print_level = 3;
    homotopy_options.sigma_0 = sigma_0;
    homotopy_options.opts_casadi_nlp.ipopt.linear_solver = 'mumps'; %
    homotopy_options.homotopy_steering_strategy = "DIRECT";
    homotopy_options.lift_complementarities = 0;

    % ccopt options
    ccopt_options.opts_madnlp.linear_solver = 'Ma27Solver';
    ccopt_options.opts_madnlp.tol = tol;
    ccopt_options.opts_madnlp.barrier.TYPE = 'MonotoneUpdate';
    ccopt_options.opts_madnlp.disable_garbage_collector = true;
    % ccopt_options.opts_madnlp.barrier.mu_min = 1e-9;
    % ccopt_options.opts_ccopt.relaxation_update.TYPE = 'RolloffRelaxationUpdate';
    % ccopt_options.opts_ccopt.relaxation_update.rolloff_slope = 2.5;
    % ccopt_options.opts_ccopt.relaxation_update.rolloff_point = 1e-12;
    % ccopt_options.opts_ccopt.sigma_min = 1e-8;
    %ccopt_options.opts_ccopt.relaxation_update.TYPE = 'RelaxLBUpdate';
    %ccopt_options.opts_ccopt.relaxation_update.relax_threshold = 1e-5;
    ccopt_options.opts_ccopt.q_regularization = 'critical_rho';
    ccopt_options.opts_ccopt.critical_rho_factor = 0.9999;
    %ccopt_options.opts_madnlp.barrier.mu_init = 1.0;
    %ccopt_options.opts_madnlp.barrier.TYPE = 'QualityFunctionUpdate';

    % MPC options
    mpc_options.fast_sigma_0 = 1; % Set the fast sigma for followup solves to the complementarity tol to only do 1 nlp solve.
    mpc_options.solve_advanced_problem = solve_advanced_problem;
    mpc_options.sqpec_hessian_convexification = "MIRROR";
    mpc_options.discard_constraints_in_hessian = 1;
    mpc_options.rho_d = 0;
    mpc_options.advanced_n_qpecs = advanced_n_qpecs;
    mpc_options.do_shift_initialization = do_shift_initialization;

    mpc_options.use_probing_qp = use_probing_qp;
    mpc_options.use_feedback_qp = use_feedback_qp; %
    mpc_options.objective_ratio = 0.97;
    mpc_options.warmstart_qpec = warmstart_qpec;
    video_name = filename_pdf;

    %% Discretization
    
    % --- fixed times (seconds) ---
    if 0
        % prototype
        DT = 0.1;
        T_pred           = 1.0;    % prediction horizon
        t_new_set_point  = 5;
        t_new_set_point  = 4;
        t_disturbance    = 4;
        T_sim            = 6.0;   % total simulated time
    else
        % expriemnts
        DT = 0.1;
        %DT = 0.1;
        T_pred           = 1.0;    % prediction horizon
        t_new_set_point  = 6;
        t_disturbance    = 10;
        T_sim            = 16.0;   % total simulated time
    end


    % --- convert to step counts (integers) ---
    N_step_new_set_point = round(t_new_set_point/DT);
    N_step_disturbance   = round(t_disturbance/DT);

    N_steps  = round(T_sim/DT);
    N_stages = round(T_pred/DT);


    %% friction parameters
    F_friction_model = 3; % Friction force amplitude
    F_friction_plant = 3;


    %
    %% nosnoc setup
    problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.

    %% Model variables and parameters
    % specify initial and desired state
    x0 = [1; 0/180*pi; 0; 0]; % start downwards
    x_ref = [0; 180/180*pi; 0; 0]; % end upwards
    x_ref1 = x_ref;
    x_ref2 = [-1; 180/180*pi; 0; 0]; % move to left
    u_ref = 0;
    % x0 = x_ref;
    model = nosnoc.model.Pss();
    m1 = 1; % cart
    m2 = 0.1; % link
    g = 9.81;
    link_length = 1;

    % CasADi symbolic variables
    n_x = 4;
    n_u = 1;
    px = SX.sym('px');
    theta = SX.sym('theta');
    v = SX.sym('v');
    theta_dot = SX.sym('theta_dot');
    q = vertcat(px, theta);
    q_dot = vertcat(v, theta_dot);
    x = vertcat(q, q_dot); % state vector
    u = SX.sym('u'); % control
    % Parameters
    F_friction = SX.sym('F_friction');
    y_ref_p = SX.sym('x_ref_p',n_x+n_u); % references

    %% Model
    M = [m1 + m2, m2*link_length*cos(theta);...
        m2 *link_length*cos(theta),  m2*link_length^2]; % Inertia matrix
    % Coriolis force
    Cor = [0, -m2 * link_length*theta_dot*sin(theta); 0,   0];
    f_all = [0; -m2*g*link_length*sin(theta)] + [u; 0] - Cor*q_dot;

    if 1
        f_0 = [q_dot; inv(M)*(f_all)];
        % Dynamics for v>0
        f_1 = [zeros(2,1);-inv(M)*[F_friction; 0]];
        % Dynamics for v<0
        f_2 = [zeros(2,1);inv(M)*[F_friction; 0]];
    else
        f_0 = zeros(4,1);
        % Dynamics for v>0
        f_1 = [q_dot;inv(M)*(f_all-[F_friction; 0])];
        % Dynamics for v<0
        f_2 = [q_dot; inv(M)*(f_all+[F_friction; 0])];
    end

    f_1_fun = Function('f_1_fun',{x,u,F_friction},{f_1});
    f_2_fun = Function('f_2_fun',{x,u,F_friction},{f_2});

    F = [f_1_fun(x,u,F_friction_model), f_2_fun(x,u,F_friction_model)];

    ubx = [10; inf; inf; inf];
    lbx = [-10; -inf; -inf; -inf];
    u_max = 20;


    % stage cost
    f_q = (x-y_ref_p(1:4))'*Q*(x-y_ref_p(1:4))+ (u-y_ref_p(5))'*R*(u-y_ref_p(5)); % running/stage costs
    % terminal cost
    f_q_T  = (x-y_ref_p(1:4))'*Q_terminal*(x-y_ref_p(1:4));

    % Populate nosnoc model
    model.c = v;         % switching function c: cart velocity
    model.S = [1; -1];   % sign matrix S % f_1 for c>0, f_2 for c<0
    model.p_global = y_ref_p;
    model.p_global_val = [x_ref;u_ref];
    model.F = F;
    model.f_0 = f_0;
    model.f_q = f_q;
    model.f_q_T = f_q_T;
    model.lbx = lbx;
    model.ubx = ubx;
    model.x = x;
    model.x0 = x0;
    model.u = u;
    model.lbu = -u_max;
    model.ubu = u_max;
    %
    % model.g_terminal = model.x-x_ref;
    % problem_options.relax_terminal_constraint = "ELL_INF";

    %% create plant
    if model_plant_missmatch
        % Integrator model
        sim_model = nosnoc.model.Pss();
        sim_model.c = v;
        sim_model.x = x;
        sim_model.x0 =  x0;
        sim_model.u = u;
        sim_model = get_cart_pole_with_friction_model(true, F_friction_plant);
        sim_model.lbx = -inf(n_x,1);
        sim_model.ubx = inf(n_x,1);
        sim_model.lbu = -inf;
        sim_model.ubu = inf;
        sim_problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
        integrator_options = nosnoc.integrator.Options();
        integrator_options.integrator_plugin = "FESD";
        %integrator_options.print_level = 0;

        sim_solver_options = integrator_options.fesd_solver_opts; % the fesd integrator uses an mpec solver, call and modify its options
        sim_problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
        if model_plant_same_integrator
            N_sim = 1;
            sim_problem_options.n_s = n_s; % Number of stage points in the RK method (determines accuracy)
            sim_problem_options.N_finite_elements = 2; % Number of finite element (integration steps) on every control interval (optionally a vector might be passed).
            sim_problem_options.N_sim = N_sim;
        else
            sim_problem_options.n_s = n_s_sim; % Number of stage points in the RK method (determines accuracy)
            sim_problem_options.N_finite_elements = 2; % Number of finite element (integration steps) on every control interval (optionally a vector might be passed).
            sim_problem_options.N_sim = N_sim;
        end
        sim_problem_options.T_sim = DT;
        sim_solver_options.print_level = 0;
        sim_solver_options.homotopy_steering_strategy = HomotopySteeringStrategy.DIRECT; % Use the $\ell_{\infty}$ steering strategy
        sim_solver_options.complementarity_tol = 1e-9; % Value to drive the complementarity residual to.

        sim_solver_options.opts_casadi_nlp.ipopt.linear_solver = 'mumps'; % Using an HSL solver is a significant speedup over the default 'mumps' but requires installation

        % create nosnoc integrator
        integrator = nosnoc.Integrator(sim_model, sim_problem_options, integrator_options);
    end
    %% OCP problem options
    problem_options.T = N_stages*DT;  % Time horizon
    problem_options.rk_scheme = RKSchemes.RADAU_IIA;
    problem_options.n_s = n_s;
    problem_options.dcs_mode = 'Stewart';
    problem_options.N_stages = N_stages; % number of control intervals
    problem_options.N_finite_elements = N_FE; % number of finite element on every control interval
    problem_options.cross_comp_mode = "FE_FE";
    problem_options.use_fesd = use_fesd;
    problem_options.gamma_h = gamma_h;
    % problem_options.lift_complementarities = 1;
    problem_options.euler_cost_integration = 1;
    % problem_options.euler_cost_integrat

    %% create mpc object
    if use_rtmpc
        if use_ccopt
            mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, ccopt_options, ccopt_options);
        else
            mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, homotopy_options, homotopy_options);
        end
        % mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, homotopy_options, mpecopt_options);
    else
        if use_ccopt
            mpc = nosnoc.mpc.FullMpcc(model, mpc_options, problem_options, ccopt_options);
        else
            mpc = nosnoc.mpc.FullMpcc(model, mpc_options, problem_options, homotopy_options);
        end
    end
    %% MPC loop
    x0 = model.x0;
    u = [];
    t = 0;
    tf = [];
    x = x0;
    f_opt = [];
    x_high_res = x0;
    t_high_res = 0;
    x_mpc_pred = x0;


    feedback_times = [];
    preparation_times = [];
    qp_solved = [];
    qp_accepted = [];
    qpec_solved = [];
    tp_i = 0; % very first prepartion;

    if plot_intermediate_solutions
        q_plot = subplot(311); hold on;
        xlim([0 (N_steps+N_stages)*problem_options.h])
        ylim([-pi*1.1 pi*1.1])
        v_plot = subplot(312); hold on;
        xlim([0 (N_steps+N_stages)*problem_options.h])
        ylim([-5 5])
        u_plot = subplot(313); hold on;
        xlim([0 (N_steps+N_stages)*problem_options.h])
        ylim([-u_max u_max])
    end

    for step=1:N_steps

        stats_prep = mpc.do_preparation();
        tp_i = stats_prep.preparation_time;
        % what would the plant to? - compute model plat missmatch for current u_mpc
        preparation_times = [preparation_times, tp_i ];

        [u_i, stats_fb] = mpc.get_feedback(x0);
        % f_opt = [f_opt, mpc.get_objective()];
        tf_i = stats_fb.feedback_time;
        feedback_times = [feedback_times, tf_i];
        if use_rtmpc
            if save_qpecs && mod(step-1, 10) == 0
                save_qpec(mpc.qpec, ['CARTPOLE_CONVEX_', num2str(step)])
            end
            qp_solved = [qp_solved, stats_fb.qp_solved];
            qp_accepted = [qp_accepted, stats_fb.qp_accepted];
            qpec_solved = [qpec_solved, stats_fb.qpec_solved];
        else
            qp_solved = [qp_solved, nan];
            qp_accepted = [qp_accepted, nan];
            qpec_solved = [qpec_solved, nan];
        end


        fprintf("MPC step: %d/%d, Preparation time: %d [s], Feedback time: %d [s]\n", step, N_steps, tp_i, tf_i);
        fprintf("Control input: %2.2f\n", u_i);

        f_opt = [f_opt, mpc.get_objective()];


        if model_plant_missmatch
            % Advance plat in time
            integrator.set_x0(x0);
            [t_grid_sim, x_sim] = integrator.simulate("u", repmat(u_i, [1,N_sim]), "x0", x0);
            x_high_res = [x_high_res,x_sim];
            t_high_res = [t_high_res, t_high_res(end) + t_grid_sim];
            x0 = x_sim(:, end);
        else
            x0 = mpc.get_predicted_state();
        end


        if step == N_step_new_set_point
            % change setpoint
            x_ref  = x_ref2;
            mpc.set_param('p_global',[],[x_ref2;u_ref]);
        end

        if step == N_step_disturbance
            % add large distrubance to angle
            x0(2) = x0(2)-1;
        end

        if plot_intermediate_solutions
            x_res = mpc.get('x');
            q_res = x_res(1:2,:);
            v_res = x_res(3:4,:);
            t_grid = mpc.get_time_grid();
            u_res = mpc.get('u');
            t_grid_u = mpc.get_control_grid();

            %
            cla(q_plot); hold on;
            cla(v_plot); hold on;
            cla(u_plot); hold on;

            plot(q_plot, t, x(1:2,:))
            plot(q_plot, t(end)+t_grid, q_res)
            yline(q_plot,x_ref(1:2),'k--')

            plot(v_plot, t, x(3:4,:))
            plot(v_plot, t(end)+t_grid, v_res)
            yline(v_plot,x_ref(3:4),'k--')

            if step>1
                stairs(u_plot, t, [u, u(:,end)]')
            end
            stairs(u_plot, t_grid_u+t(end), [u_res, u_res(:,end)]')
            yline(u_plot,u_ref,'k--')

        end
        % if u_i < -19.5
        %     keyboard
        % end
        x = [x, x0];
        u = [u,u_i];
        t = [t, t(end) + problem_options.h];
    end

    %% Plot
    if model_plant_missmatch
        t_u = t;
        t = t_high_res;
        x = x_high_res;
    else
        % copy for convienece in result evaluations
        t_u = t;
        t_high_res = t;
        x_high_res = x;

    end
    %%
    f = figure;
    % detect zero crossings of x(3,:)
    v = x(3,:);
    eps_v = 1e-1;                 % <-- tune (e.g. 1e-4 ... 1e-2 depending on scale)
    state = zeros(size(v));       % -1 (neg), 0 (deadband), +1 (pos)
    state(v >  eps_v) =  1;
    state(v < -eps_v) = -1;
    % compress: ignore deadband by carrying last nonzero state forward
    s = state;
    last = 0;
    for i = 1:numel(s)
        if s(i) ~= 0
            last = s(i);
        else
            s(i) = last;
        end
    end
    % detect sign changes of the filtered sign
    idx_cross = find(s(1:end-1).*s(2:end) < 0);
    t_cross = t(idx_cross);

    subplot(311)
    plot(t,x(1,:),'LineWidth',linewidth,'DisplayName','$q$')
    hold on
    plot(t,x(2,:),'LineWidth',linewidth,'DisplayName','$\theta$')
    % xline(DT*N_step_disturbance,'-','Color',matlab_magenta,'LineWidth',linewidth)
    % xline(DT*N_step_new_set_point,':','Color',matlab_orange,'LineWidth',linewidth)
    xlabel("$t$")
    ylabel("$q,\, \theta$")
    % grid on
    t_switch_set_point = DT*N_step_new_set_point;
    i_sw = find(t <= t_switch_set_point, 1, 'last');
    t1 = t(1:i_sw);
    t2 = t(i_sw:end);
    % reference, first segment
    plot(t1, x_ref1(1)*ones(size(t1)),'k--','LineWidth',linewidth,'HandleVisibility','off')
    plot(t1, x_ref1(2)*ones(size(t1)),'k--','LineWidth',linewidth,'HandleVisibility','off')

    % reference, second segment
    plot(t2, x_ref2(1)*ones(size(t2)),'k--','LineWidth',linewidth,'HandleVisibility','off')
    plot(t2, x_ref2(2)*ones(size(t2)),'k--','LineWidth',linewidth,'HandleVisibility','off')

    legend('Location','southeast');
    %
    % yline(x_ref1(1:2),'k--','LineWidth',1.5)
    % yline(x_ref2(1:2),'k-','LineWidth',1.5)

    subplot(312)
    plot(t,x(3,:),'LineWidth',linewidth,'DisplayName','$v$')
    hold on
    plot(t,x(4,:),'LineWidth',linewidth,'DisplayName','$\omega$')
    for k = 1:numel(t_cross)
        xline(t_cross(k),'k:','LineWidth',1.2,'HandleVisibility','off')
    end
    xlabel("$t$")
    ylabel("$v,\, \omega$")
    legend('Location','southeast');
    % grid on

    subplot(313)
    stairs(t_u,[u,u(end)],'LineWidth',linewidth)
    xlabel("$t$")
    ylabel("$u$")
    grid on
    ylim([-1.1*u_max 1.1*u_max])
    hold on
    yline(-u_max,'--','Color',matlab_blood_red,'HandleVisibility','off')
    yline(u_max,'--','Color',matlab_blood_red,'HandleVisibility','off')
    % exportgraphics(f, [ filename_pdf '1.pdf']);



    %%
    % figure
    % subplot(121)
    % semilogy(t,vecnorm(x-x_ref));
    % xlabel("$t$")
    % ylabel("$|x-x_{\mathrm{sp}}|$")
    % grid on
    % subplot(122)
    % semilogy(t,[nan,f_opt])
    % xlabel("$t$")
    % ylabel("$f^*$ (objective of every mpc subproblem)")
    % grid on

    %% Tracking error plots for single trajectory (x - x_sp(t))
    f = figure;

    t_switch_set_point = t_new_set_point;

    % build piecewise reference over t
    xref = repmat(x_ref1, 1, numel(t));
    idx2 = t > t_switch_set_point;
    xref(:,idx2) = repmat(x_ref2, 1, nnz(idx2));

    e = abs(x - xref);   % tracking error

    subplot(411)
    semilogy(t,e(1,:),'LineWidth',linewidth,'DisplayName','$e_q$'); hold on
    semilogy(t,e(2,:),'LineWidth',linewidth,'DisplayName','$e_\theta$');

    xlabel("$t$")
    ylabel('$\|q-q_{sp}\|_2$')
    yline(0,'k--','HandleVisibility','off')
    legend('Location','southeast');

    subplot(412)
    semilogy(t,e(3,:),'LineWidth',linewidth,'DisplayName','$e_v$'); hold on
    semilogy(t,e(4,:),'LineWidth',linewidth,'DisplayName','$e_\omega$');
    xlabel("$t$")
    ylabel('$\|v-v_{sp}\|_2$')
    yline(0,'k--','HandleVisibility','off')
    legend('Location','southeast');

    subplot(413)
    err_norm = vecnorm(e,2,1);
    semilogy(t, err_norm,'LineWidth',linewidth,'DisplayName','$\|x-x_{sp}\|_2$');
    xlabel("$t$")
    ylabel('$\|x-x_{sp}\|_2$')
    grid on
    legend('Location','best');

    subplot(414)
    semilogy(f_opt,'LineWidth',linewidth);
    xlabel("$k$")
    ylabel('$V(x^k)$')
    grid on






    %% Animation
    N_history = round(0.1*length(x));
    if create_video
        outputVideo = VideoWriter([video_name '_cart.mp4'], 'MPEG-4');
        outputVideo.FrameRate = 25; % Adjust the frame rate as needed
        open(outputVideo);
        time_step = mean(diff(t));
        q1_opt = x(1,:);
        q2_opt = x(2,:);
        v1_opt = x(3,:);
        v2_opt = x(4,:);
        t_grid = t;
        t_grid_u = t;
        u_opt = [u,u(end)];
        % time_step = problem_options.h;
        % filename = 'cart_pole_with_friction.gif';
        figure('Renderer', 'painters', 'Position', [100 100 1200 600])
        % figure
        link_length = 1;
        cart_center = 0.125;
        cart_width1 = 0.25;
        cart_height = cart_center*2;
        pole_X = [q1_opt',q1_opt'+(link_length)*cos(q2_opt'-pi/2)];
        pole_Y = [cart_center+0*q1_opt',cart_center+link_length*sin(q2_opt'-pi/2)];
        x_min =-3;
        x_max = 3;
        for ii = 1:length(q1_opt)
            % pole
            plot(pole_X(ii,:),pole_Y(ii,:),'k','LineWidth',3);
            hold on
            % tail
            history_index = max(1,ii-N_history):ii;
            plot(pole_X(history_index,2),pole_Y(history_index,2),'color',[1 0 0 0.5],'LineWidth',0.5);
            % cart
            xp = [q1_opt(ii)-cart_width1/2 q1_opt(ii)+cart_height/2 q1_opt(ii)+cart_height/2 q1_opt(ii)-cart_width1/2];
            yp = [0 0 cart_height  cart_height];
            patch(xp,yp,'k','FaceAlpha',0.8)

            % targent 1
            % pole
            plot([x_ref(1),x_ref(1)+(link_length)*cos(x_ref(2)-pi/2)],...
                [cart_center+0,cart_center+link_length*sin(x_ref(2)-pi/2)],'color',[0 0 0 0.1],'LineWidth',3);
            % cart
            xp = [x_ref(1)-cart_width1/2 x_ref(1)+cart_height/2 x_ref(1)+cart_height/2 x_ref(1)-cart_width1/2];
            yp = [0 0 cart_height  cart_height];
            patch(xp,yp,'k','FaceAlpha',0.1)

            % targent 1
            % pole
            plot([x_ref1(1),x_ref1(1)+(link_length)*cos(x_ref1(2)-pi/2)],...
                [cart_center+0,cart_center+link_length*sin(x_ref1(2)-pi/2)],'color',[0 0 0 0.15],'LineWidth',3);
            % cart
            xp = [x_ref1(1)-cart_width1/2 x_ref1(1)+cart_height/2 x_ref1(1)+cart_height/2 x_ref1(1)-cart_width1/2];
            yp = [0 0 cart_height  cart_height];

            patch(xp,yp,'k','FaceAlpha',0.15)

            x_min = -3+q1_opt(ii);
            x_max = +3+q1_opt(ii);
            % ground
            xp = [x_min x_max x_max x_min ];
            yp = [-2 -2 0 0];
            patch(xp,yp,0*ones(1,3),'FaceAlpha',0.1,'EdgeColor','none');

            axis equal
            xlim([x_min x_max])
            ylim([-1 2])
            text(-1.5,1.5,['Time: ' num2str(t_grid(ii),'%.2f') ' s'],'interpreter','latex','fontsize',15)

            frame = getframe(gcf);
            writeVideo(outputVideo, frame);
            clf;
        end
        close(outputVideo);
    end

    %%
    %% SAVE DATA (single struct, file name = data_name.mat)
    if save_data
        data = struct();
        data.name = data_name;
        data.scenario = scenario;

        % trajectories
        data.t = t;               % state time grid (high-res if mismatch)
        data.x = x;
        data.t_u = t_u;           % control time grid (MPC grid)
        data.u = u;

        data.f_opt = f_opt; % value function

        % stats
        data.preparation_times = preparation_times;
        data.feedback_times = feedback_times;

        % high-res explicitly (always stored)
        data.t_high_res = t_high_res;
        data.x_high_res = x_high_res;

        % timing / solver stats
        data.tf = tf;

        % references + event times
        data.x_ref1 = x_ref1;
        data.x_ref2 = x_ref2;
        data.u_ref = u_ref;
        data.t_new_set_point = t_new_set_point;
        data.t_disturbance = t_disturbance;

        % settings (useful for bookkeeping)
        data.DT = DT;
        data.T_pred = T_pred;
        data.T_sim = T_sim;
        data.N_steps = N_steps;
        data.N_stages = N_stages;
        data.N_sim = N_sim;
        data.u_max = u_max;
        data.F_friction_model = F_friction_model;
        data.F_friction_plant = F_friction_plant;
        data.model_plant_missmatch = model_plant_missmatch;
        data.use_rtmpc = use_rtmpc;
        data.gamma_h = gamma_h;
        data.N_FE = N_FE;
        data.n_s = n_s;
        data.n_s_sim = n_s_sim;

        save([data_name '.mat'], 'data');
        fprintf('Saved dataset to %s.mat\n', data_name);
    end
end
