%% spinner_mpc.m 
latexify_plot()
nice_plot_colors;
clear; clc; close all;
import casadi.*

%% setup
filename_pdf = 'cartpole_mpc';
scenarios= {'Full-MPCC-MPC','HyRTI','HyAS-RTI-QPCC','HyAS-RTI-QPCC-FULL','HyRTI-IPOPT','HyAS-RTI-QPCC-IPOPT','Smoothing','HyRTI-fast','HyAS-RTI-QP'};
%% ---------------- Settings ----------------
shape = 'circle';                   % 'circle' or 'ellipse'
run_animation = 1;
filename = 'spinner_mpc_traj';
use_rtmpc = 1;
use_ccopt = 1;
save_qpecs = 0;
save_data = 1;
model_plant_mismatch = 1;
n_s_sim = 1;
N_sim = 8;
for kk = [7]%[1,2,3,4,7]%1:length(scenarios)
    scenario =  scenarios{kk};
    % initalize opts
    mpc_options = nosnoc.mpc.Options(); % Initialize all options related to the MPC implementation.
    homotopy_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosonc problems.
    ccopt_options = nosnoc.ccopt.Options(); % ccopt options
    ccopt_options_qpcc = nosnoc.ccopt.Options(); % ccopt options
    problem_options = nosnoc.Options();
    % MPEC solver options
    tol = 1e-6;
    qpcc_mu_init = 1e-1;
    advanced_problem_type = 'full';
    N_homotopy = 10;
    sigma_0 = 1;
    fast_sigma_0 = 1e-3; % sigma for mpec warmstart in mpc

    switch scenario
        case 'Full-MPCC-MPC'
            data_name = scenario;
            use_rtmpc = 0;
            do_shift_initialization = 0;
            solve_advanced_problem = 0;
            advanced_n_qpecs = 1;
            use_feedback_qp = false;
            use_probing_qp = false;
            warmstart_qpec = false;
            tol = 1e-6;
        case 'HyRTI'
            data_name = scenario;
            use_rtmpc = 1;
            use_ccopt = 1;
            do_shift_initialization = 0;
            solve_advanced_problem = 0;
            advanced_n_qpecs = 1;
            use_feedback_qp = false;
            use_probing_qp = false;
            warmstart_qpec = true;
            qpcc_mu_init = 1e-1;
            % N_homotopy = 5;
            % sigma_0 = 1;
            tol = 1e-6;
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
            tol = 1e-6;
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
        case 'HyAS-RTI-QPCC-FULL'
            data_name = scenario;
            use_rtmpc = 1;
            use_ccopt = true;
            do_shift_initialization = 0;
            solve_advanced_problem = 1;
            advanced_problem_type = 'full';
            advanced_n_qpecs = 3;
            use_feedback_qp = false;
            use_probing_qp = false;
            warmstart_qpec = true;
            qpcc_mu_init = 1e-1;
            tol=1e-6;
        case 'HyAS-RTI-QPCC'
            data_name = scenario;
            use_rtmpc = 1;
            use_ccopt = true;
            do_shift_initialization = 0;
            solve_advanced_problem = 1;
            advanced_problem_type = 'sqpec';
            advanced_n_qpecs = 1;
            use_feedback_qp = false;
            use_probing_qp = false;
            warmstart_qpec = true;
            qpcc_mu_init = 1e-1;
            tol=1e-6;
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
            use_ccopt = false;
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
    %% Discretization

    % --- fixed times (seconds) ---
    if 0
        % prototype
        DT = 0.05;
        T_pred           = .75;    % prediction horizon
        t_new_set_point  = 5;
        t_new_set_point  = 4;
        t_disturbance    = 4;
        T_sim            = 16.0;   % total simulated time
    else
        % expriemnts
        DT = 0.1;
        %DT = 0.1;
        T_pred           = 1.5;    % prediction horizon
        t_disturbance    = 6;
        t_disturbance_2  = 11;
        T_sim            = 16.0;   % total simulated time
    end


    % --- convert to step counts (integers) ---

    N_steps  = round(T_sim/DT);
    N_stages = round(T_pred/DT);

    %% --------------- NOSNOC options --------------
    problem_options.N_stages = N_stages;
    problem_options.N_finite_elements = ones(N_stages,1);
    problem_options.N_finite_elements(1) = 4;
    problem_options.N_finite_elements(2) = 4;
    problem_options.N_finite_elements(3) = 2;
    problem_options.N_finite_elements(4) = 2;
    problem_options.N_finite_elements(5) = 1;
    
    problem_options.h = DT;
    problem_options.T = T_pred;

    problem_options.rk_scheme = RKSchemes.RADAU_IIA;
    problem_options.n_s = 1;
    problem_options.cross_comp_mode = 7;
    problem_options.use_fesd = 0;
    problem_options.lift_velocity_state = 1;
    problem_options.friction_model = 'Polyhedral';
    problem_options.euler_cost_integration = false;
    problem_options.g_path_at_fe = true;

    problem_options.ub_gamma_d = 1e3;
    problem_options.ub_delta_d = 1e3;

    homotopy_options.N_homotopy = N_homotopy;
    homotopy_options.homotopy_update_slope = 0.1;
    homotopy_options.complementarity_tol = 1e-6;
    homotopy_options.sigma_0 = sigma_0;
    % solver_options.homotopy_update_rule = 'superlinear';
    homotopy_options.relaxation_strategy = "SCHOLTES_INEQ";
    homotopy_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
    homotopy_options.opts_casadi_nlp.ipopt.hsllib = '/home/anton/tools/HSL_jll.jl-2023.11.7/override/lib/x86_64-linux-gnu-libgfortran5/libhsl.so';
    homotopy_options.opts_casadi_nlp.ipopt.max_iter = 3000;
    homotopy_options.print_level = 3;

    ccopt_options.opts_madnlp.linear_solver = 'Ma27Solver';
    ccopt_options.opts_madnlp.tol = tol;
    ccopt_options.opts_madnlp.barrier.TYPE = 'MonotoneUpdate';
    ccopt_options.opts_madnlp.barrier.mu_min = 1e-1*tol;
    ccopt_options.opts_madnlp.disable_garbage_collector = true;
    ccopt_options.opts_madnlp.print_level = 6;
    % ccopt_options.opts_ccopt.relaxation_update.TYPE = 'RolloffRelaxationUpdate';
    % ccopt_options.opts_ccopt.relaxation_update.rolloff_point = 1e-3;
    % ccopt_options.opts_ccopt.relaxation_update.rolloff_slope = 1.5;
    % ccopt_options.opts_ccopt.relaxation_update.sigma_min = tol;
    ccopt_options.opts_ccopt.relaxation_update.TYPE = 'ProportionalRelaxationUpdate';
    ccopt_options.opts_ccopt.relaxation_update.sigma_min = tol;

    ccopt_options_qpcc.opts_madnlp.linear_solver = 'Ma27Solver';
    ccopt_options_qpcc.opts_madnlp.tol = tol;
    ccopt_options_qpcc.opts_madnlp.barrier.TYPE = 'MonotoneUpdate';
    ccopt_options_qpcc.opts_madnlp.barrier.mu_init = qpcc_mu_init;
    ccopt_options_qpcc.opts_madnlp.barrier.mu_min = 1e-1*tol;
    ccopt_options_qpcc.opts_madnlp.disable_garbage_collector = true;
    ccopt_options_qpcc.opts_madnlp.print_level = 6;
    % ccopt_options_qpcc.opts_ccopt.relaxation_update.TYPE = 'RolloffRelaxationUpdate';
    % ccopt_options_qpcc.opts_ccopt.relaxation_update.rolloff_point = 1e-4;
    % ccopt_options_qpcc.opts_ccopt.relaxation_update.rolloff_slope = 1.5;
    % ccopt_options_qpcc.opts_ccopt.relaxation_update.sigma_min = tol;
    ccopt_options_qpcc.opts_ccopt.relaxation_update.TYPE = 'ProportionalRelaxationUpdate';
    ccopt_options_qpcc.opts_ccopt.relaxation_update.sigma_min = tol;
    ccopt_options_qpcc.opts_ccopt.q_regularization = 'critical_rho';
    ccopt_options_qpcc.opts_ccopt.critical_rho_factor = 0.999;
    ccopt_options_qpcc.opts_ccopt.min_reg_mu = 0.0;
    %ccopt_options_qpcc.opts_ccopt.print_level = 6;

    %% MPC options
    mpc_options.fast_sigma_0 = 1; % Set the fast sigma for followup solves to the complementarity tol to only do 1 nlp solve.
    mpc_options.solve_advanced_problem = solve_advanced_problem;
    %mpc_options.sqpec_hessian_convexification = "PROJECT";
    mpc_options.discard_constraints_in_hessian = 1;
    mpc_options.rho_d = 0;
    mpc_options.advanced_problem_type = advanced_problem_type;
    mpc_options.advanced_n_qpecs = advanced_n_qpecs;
    mpc_options.do_shift_initialization = do_shift_initialization;
    mpc_options.warmstart_full_mpc = true;

    mpc_options.use_probing_qp = use_probing_qp;
    mpc_options.use_feedback_qp = use_feedback_qp; %
    mpc_options.objective_ratio = 0.97;
    mpc_options.warmstart_qpec = warmstart_qpec;
    video_name = filename_pdf;
    %% --------------- Model -----------------------
    [model, robot_data] = spinner_model(shape);
    n_x = length(model.x);
    n_q = robot_data.n_q; n_u = robot_data.n_u;
    x = model.x; u = model.u;
    model.D_tangent = [model.J_tangent, -model.J_tangent];

    q_init = [0*pi/32; 3*pi/8; 0.0];
    v_init = [0;0;0];
    x0 = [q_init; v_init];

    % Parametric reference (p_global)
    y_ref_p = SX.sym('y_ref_p', length(x)+n_u+2);

    % Weights
    Qq = diag([0.5, 0.1, 20.0]);
    Qv = diag([1.0 1.0 0.001]);
    Qcart = diag([0.0 0.0]);
    Q  = blkdiag(Qq, Qv);
    R  = diag([0.1 0.1])*1e-1;
    Qfq = diag([1 1 10.0]);
    Qfv = diag([0.1 0.1 1.0]);
    Q_terminal = blkdiag(Qfq, Qfv);

    x_ref_sym = SX.sym('x_ref',n_x);
    x_ref_end_sym = y_ref_p(1:length(x));
    cart_ref_sym = y_ref_p(length(x)+(1:2));
    u_ref_sym = y_ref_p(length(x)+2+(1:n_u));

    % Costs
    f_q = (x - x_ref_sym).' * Q * (x - x_ref_sym) + (u - u_ref_sym).' * R * (u - u_ref_sym) + (robot_data.p2 - cart_ref_sym).' * Qcart * (robot_data.p2 - cart_ref_sym);
    f_q_T = (x - x_ref_sym).' * Q_terminal * (x - x_ref_sym);

    %g_path = [robot_data.p2;x(1:2); x(4:5)];
    %lbg_path = [0.7;0.5;-pi;-pi;-30;-30];
    %ubg_path = [2.0;2.0;pi;pi;30;30];
    % g_path = [robot_data.p2];
    % lbg_path = [0.0;0.0];
    % ubg_path = [1.4;2.5];
    g_path = [robot_data.p2(1)];
    lbg_path = [-inf];
    ubg_path = [1.4];
    model.g_path = g_path;
    model.lbg_path = lbg_path;
    model.ubg_path = ubg_path;
    problem_options.relax_path_constraints = "ELL_2";
    problem_options.rho_path = 0.1;

    % Default reference
    q_nom_start = [q_init(1:2); 0];
    v_nom_zero  = [0;0;0];
    cart_target = [1.2;1.0];
    x_ref_default = [q_nom_start; v_nom_zero; cart_target];
    u_ref_default = zeros(n_u,1);
    u1_max = 15;
    u2_max = 15;
    model.ubu = [u1_max;u2_max];
    model.lbu = [-u1_max;-u2_max];
    model.x0 = x0;
    model.f_q = f_q;
    model.f_q_T = f_q_T;
    model.p_global = [x_ref_end_sym; u_ref_sym];
    model.p_time_var = [x_ref_sym;cart_ref_sym];
    model.p_global_val = [x_ref_default(1:n_x); u_ref_default];
    model.p_time_var_val = repmat(x_ref_default, 1, N_stages);
    
    %% create plant
    if model_plant_mismatch
        [sim_model, sim_robot_data] = spinner_model(shape);
        % Integrator model
        sim_model.lbx = -inf(n_x,1);
        sim_model.ubx = inf(n_x,1);
        sim_model.lbu = -inf(n_u,1);
        sim_model.ubu = inf(n_u,1);
        sim_model.D_tangent = [sim_model.J_tangent, -sim_model.J_tangent];
        sim_problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
        integrator_options = nosnoc.integrator.Options();
        integrator_options.integrator_plugin = "FESD";
        %integrator_options.print_level = 0;

        sim_solver_options = integrator_options.fesd_solver_opts; % the fesd integrator uses an mpec solver, call and modify its options
        sim_problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme

        sim_problem_options.n_s = 1; % Number of stage points in the RK method (determines accuracy)
        sim_problem_options.N_finite_elements = 1; % Number of finite element (integration steps) on every control interval (optionally a vector might be passed).
        sim_problem_options.N_sim = N_sim;
        sim_problem_options.T_sim = DT;
        sim_problem_options.use_fesd = false;
        sim_problem_options.ub_gamma_d = 1e2;
        sim_problem_options.ub_delta_d = 1e2;
        sim_problem_options.lift_velocity_state = 1;
        sim_problem_options.friction_model = 'Polyhedral';


        sim_solver_options.print_level = 1;
        sim_solver_options.homotopy_steering_strategy = HomotopySteeringStrategy.ELL_INF; % Use the $\ell_{\infty}$ steering strategy
        sim_solver_options.complementarity_tol = 1e-7; % Value to drive the complementarity residual to.

        %sim_solver_options.opts_casadi_nlp.ipopt.linear_solver = 'mumps'; % Using an HSL solver is a significant speedup over the default 'mumps' but requires installation
        sim_solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';
        sim_solver_options.opts_casadi_nlp.ipopt.hsllib = '/home/anton/tools/HSL_jll.jl-2023.11.7/override/lib/x86_64-linux-gnu-libgfortran5/libhsl.so';
        % create nosnoc integrator
        integrator = nosnoc.Integrator(sim_model, sim_problem_options, integrator_options);
    end
    %% --------------- Build MPC solver -------------
    if use_rtmpc
        if use_ccopt
            mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, ccopt_options, ccopt_options_qpcc);
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

    %% --------------- Reference schedule -----------
    % Edit these two vectors only
    theta_ref_seq      = [-1*pi/4; pi/2; -pi/2; 1*pi/4; 0];    % target angles (rad)
    theta_change_steps = [1, 30, 80, 110. 140];        % MPC steps where they become active (1-based)
    
    theta_ref_seq      = [0; -2*pi/4; -2*pi/4; 3*pi/8; 3*pi/8; -pi/2; -pi/2; pi/4; pi/4;   0;   0;   0];    % target angles (rad)
    theta_change_steps = [1,      20,      30,     50,     70,    95,   120,  135,  145, 146, 160, 176];        % MPC steps where they become active (1-based)

    % Safety: ensure both start at step 1
    if theta_change_steps(1) ~= 1
        theta_change_steps = [1, theta_change_steps];
        theta_ref_seq      = [theta_ref_seq(1); theta_ref_seq];
    end
    assert(numel(theta_ref_seq) == numel(theta_change_steps), 'Lengths must match.');

    % Build full state references per change
    x_refs = repmat([q_nom_start(1:2); 0; v_nom_zero;cart_target], 1, N_steps+N_stages+1);
    interp_ref = 'linear';
    for i = 1:(N_steps+N_stages+1)
        x_refs(3,i) = interp1(theta_change_steps, theta_ref_seq, i, interp_ref);
    end

    % Initialize
    mpc.set_param('p_global', [], [x_refs(1:n_x,N_stages+1); u_ref_default]);
    for stage = 1:N_stages
        mpc.set_param('p_time_var', {stage}, [x_refs(:,stage)]);
    end
    current_ref_idx = 1;

    %% --------------- Closed-loop simulation -------
    x0 = model.x0;
    u = [];
    t = 0;
    tf = [];
    x = x0;
    f_opt = [];
    x_high_res = x0;
    t_high_res = 0;
    x_mpc_pred = x0;
    lambda_n_high_res = [0];


    feedback_times = [];
    preparation_times = [];
    qp_solved = [];
    qp_accepted = [];
    qpec_solved = [];
    for step=1:N_steps
        stats_prep = mpc.do_preparation();
        tp_i = stats_prep.preparation_time;
        % what would the plant to? - compute model plat missmatch for current u_mpc
        preparation_times = [preparation_times, tp_i ];

        % Online reference update
        mpc.set_param('p_global', [], [x_refs(1:n_x,N_stages+step+1); u_ref_default]);
        for stage = 1:N_stages
            mpc.set_param('p_time_var', {stage}, [x_refs(:,step+stage)]);
        end
        
        [u_i, stats_fb] = mpc.get_feedback(x0);
        tf_i = stats_fb.feedback_time;
        feedback_times = [feedback_times, tf_i];
        if use_rtmpc
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


        if model_plant_mismatch
            % Advance plat in time
            integrator.set_x0(x0);
            [t_grid_sim, x_sim] = integrator.simulate("u", repmat(u_i, [1,N_sim]), "x0", x0);
            x_high_res = [x_high_res,x_sim];
            t_high_res = [t_high_res, t_high_res(end) + t_grid_sim];
            lambda_n_high_res = [lambda_n_high_res, lambda_n_high_res(end), integrator.get("lambda_normal")];
            x0 = x_sim(:, end);
        else
            x0 = mpc.get_predicted_state();
        end

        t = [t, t(end) + problem_options.h];
        % if abs(t(end) - t_disturbance) <= 1e-5
        %     x0(6) = -5;
        % end
        if abs(t(end) - t_disturbance_2) <= 1e-5
            x0(3) = x0(3) + pi/4;
        end
        x = [x, x0];
        u = [u,u_i];
        norm(x0 - mpc.get_predicted_state(),2)
        %keyboard;
    end

    %% Extract Data
    if model_plant_mismatch
        t_u = t;
        t = t_high_res;
        x = x_high_res;
    else
        % copy for convienece in result evaluations
        t_u = t;
        t_high_res = t;
        x_high_res = x;
    end
    %% Plot
    q_res = x(1:n_q,:);   % size: n_q x (N_steps+1)
    v_res = x(n_q+1:end,:);
    % --------------- Plot tracking (fixed stairs) ----------------
    figure('Color','w');

    % Actual theta over full time vector (align lengths)
    plot(t, q_res(3,:), 'LineWidth', 1.8); hold on;

    % Build a proper piecewise-constant reference over [0, t_end]
    t_end = t(end);
    t_changes = DT*(theta_change_steps - 1);   % time of each change
    % Ensure coverage starts at 0 exactly
    if t_changes(1) ~= 0
        t_changes = [0, t_changes];
        theta_ref_seq = [theta_ref_seq(1); theta_ref_seq];
    end
    % Append final time so the last step extends to t_end
    t_ref_plot = [t_changes, t_end];
    theta_ref_plot = [theta_ref_seq.', theta_ref_seq(end)];
    
    stairs(t_u(1:end-1), x_refs(3,1:(length(t_u)-1)), 'r--', 'LineWidth', 1.6);

    xlabel('t [s]'); ylabel('\theta (rad)');
    legend({'Actual \theta','Reference'}, 'Location','Best');
    title('Spinner angle tracking');
    grid on;

    %% --------------- Optional animation ----------
    if run_animation
        h = figure('Color','w'); hold on; axis equal
        xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
        xlim([-0.2 2.2]); ylim([-0.2 2.0]);

        l1 = robot_data.l1; l2 = robot_data.l2;
        cx = robot_data.cxy(1); cz = robot_data.cxy(2);
        N_frames = size(q_res,2);

        %outputVideo = VideoWriter([filename '.mp4'],'MPEG-4');
        outputVideo = VideoWriter([filename '.avi']);
        outputVideo.FrameRate = 1/(DT/N_sim); open(outputVideo);

        for k = 1:N_frames
            cla
            q1 = q_res(1,k); q2 = q_res(2,k); q3 = q_res(3,k);
            p0 = [0;0];
            p1 = p0 + [l1*sin(q1); l1*cos(q1)];
            p2 = p1 + [l2*sin(q1+q2); l2*cos(q1+q2)];
            p_ball = robot_data.cxy;

            switch lower(robot_data.shape)
                case 'circle'
                    R = robot_data.a;
                    nSlices = 4;
                    colors = [0.2 0.6 1.0;   % blue
                        1.0 0.4 0.3;   % red
                        0.2 0.6 1.0;   % blue again
                        1.0 0.4 0.3];  % red again (diagonal match)

                    % each slice is 90°
                    for s = 1:nSlices
                        th1 = (s-1)*pi/2;
                        th2 = s*pi/2;
                        th  = linspace(th1,th2,50);
                        % rotate entire circle by q3
                        rot = [cos(q3) -sin(q3); sin(q3) cos(q3)];
                        pts = rot * [ [0, R*cos(th), 0]; [0, R*sin(th), 0] ];
                        fill(cx + pts(1,:), cz + pts(2,:), colors(s,:), ...
                            'EdgeColor','none','FaceAlpha',0.8);
                    end

                    % black outline and rotation marker
                    th = linspace(0,2*pi,100);
                    plot(cx + R*cos(th), cz + R*sin(th), 'k','LineWidth',1.5);
                    plot(cx, cz, 'ko','MarkerFaceColor','k','MarkerSize',5);
                    % add a small pointer
                    % plot([cx, cx + 0.9*R*cos(q3)], [cz, cz + 0.9*R*sin(q3)], ...
                    %     'k-', 'LineWidth', 2);

                case 'ellipse'
                    a = robot_data.a; b = robot_data.b;
                    th = linspace(0, 2*pi, 200);
                    R2 = [cos(q3) -sin(q3); sin(q3) cos(q3)];
                    P = R2 * [a*cos(th); b*sin(th)];
                    plot(P(1,:)+cx, P(2,:)+cz, 'b', 'LineWidth', 2);
            end

            plot([p0(1) p1(1)], [p0(2) p1(2)], 'r','LineWidth',6);
            plot([p1(1) p2(1)], [p1(2) p2(2)], 'r','LineWidth',6);
            plot(p0(1), p0(2), 'ko','MarkerFaceColor','k');
            plot(p1(1), p1(2), 'ko','MarkerFaceColor','k');
            plot(p2(1), p2(2), 'ko','MarkerFaceColor','k');
            plot(p_ball(1), p_ball(2), 'bo','MarkerFaceColor','b');

            title(sprintf('t = %.2f s', t_high_res(k)));
            drawnow
            writeVideo(outputVideo, getframe(h));
        end
        close(outputVideo);
    end
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
        data.lambda_n_high_res = lambda_n_high_res;

        % Ref
        data.x_refs = x_refs(:,1:N_steps);

        % timing / solver stats
        data.tf = tf;

        % settings (useful for bookkeeping)
        data.DT = DT;
        data.T_pred = T_pred;
        data.T_sim = T_sim;
        data.N_steps = N_steps;
        data.N_stages = N_stages;
        data.N_sim = N_sim;
        data.u_max = [u1_max;u2_max];
        data.model_plant_mismatch = model_plant_mismatch;
        data.use_rtmpc = use_rtmpc;
        data.N_FE = problem_options.N_finite_elements;
        data.n_s = 1;
        data.n_s_sim = n_s_sim;

        save([data_name '.mat'], 'data');
        fprintf('Saved dataset to %s.mat\n', data_name);
    end
end

%% qp probing
% check:
N_check = N_steps-1;
fprintf(['QP solved: %.1f%%  |  QP accepted: %.1f%%  |  QPEC solved (fallback): %.1f%%\n'], 100*sum(qp_solved)/N_check, 100*sum(qp_accepted)/N_check, 100*sum(qpec_solved)/N_check);


%% todo
% plot controls and contact forces;
