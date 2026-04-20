% One manipulator disc pushes one passive disc to a target.
% Based on PDSObjects (NOSNOC)

clear all; close all;
import casadi.*

%% Create model and options
model = nosnoc.model.PDSObjects();
problem_options = nosnoc.Options();
reg_solver_options  = nosnoc.reg_homotopy.Options();
ccopt_options  = nosnoc.ccopt.Options();
ccopt_options_2  = nosnoc.ccopt.Options();
mpc_options = nosnoc.mpc.Options(); % Initialize all options related to the MPC implementation.

%% Parameters
T = 4.0;         % horizon (shorter, simpler)
R = 1;           % manipulator radius
R_obj = 2;       % object radius
n_d = 2;
N_steps = 100;
N_stages = 30;
N_fe = 2;
use_rtmpc = 0; % Requires CCOpt initialization! 
model_plant_mismatch = 0;

%% mpecopt as qpcc
mpecopt_options = mpecopt.Options();
mpecopt_options.compute_tnlp_stationary_point = false;
mpecopt_options.settings_lpec.lpec_solver = "Gurobi";
mpecopt_options.settings_casadi_nlp.ipopt.linear_solver = 'ma27';
mpecopt_options.relax_and_project_homotopy_parameter_steering = "Direct";
mpecopt_options.rho_TR_phase_i_init = 1e-4;
mpecopt_options.rho_TR_phase_ii_init = 1e-6;
% mpecopt_options.allow_early_termination = true;
mpc_options.objective_ratio = 0.96;

ccopt_options_2.opts_madnlp.linear_solver = 'Ma27Solver';
ccopt_options_2.opts_madnlp.tol = 1e-6;
ccopt_options_2.opts_madnlp.barrier.TYPE = 'MonotoneUpdate';
ccopt_options_2.opts_madnlp.barrier.mu_min = 1e-9;
ccopt_options_2.opts_ccopt.relaxation_update.TYPE = 'RolloffRelaxationUpdate';
ccopt_options_2.opts_ccopt.relaxation_update.rolloff_slope = 2.5;
ccopt_options_2.opts_ccopt.relaxation_update.rolloff_point = 1e-12;

ccopt_options.opts_madnlp.linear_solver = 'Ma27Solver';
ccopt_options.opts_madnlp.tol = 1e-6;
% ccopt_options.opts_ccopt.sigma_min = 1e-8;

%% Controls
u1 = SX.sym('u1', n_d);
model.u  = u1;
model.u0 = [0;0];
model.lbu = [-3;-3];
model.ubu = [ 3; 3];

y_ref_p = SX.sym('y_ref_p',n_d*2); % references


%% Objects
manip = nosnoc.objects.Ball(R, n_d);      % actuated
obj   = nosnoc.objects.Ball(R_obj, n_d);  % passive

% Initial states
manip.x0 = [-5; 5];
obj.x0   = [ 0; 0];

% Dynamics
manip.f_rhs = u1;        % actuated
obj.f_rhs   = [0; 0];    % passive

% (Optional) simple bounds to keep things in view
manip.lbx = [-20; -20];  manip.ubx = [20; 20];
obj.lbx   = [-20; -20];  obj.ubx   = [20; 20];

%% Contact: manipulator ↔ object
model.addContact(manip, obj);

%% Target & costs (focus on moving the object to y=8)
manip_target = [-5; 5];        % keep manip near start (weakly)
obj_target   = [ 0; 4];        % goal for the passive object
x_target = [manip_target; obj_target];

model.p_global = y_ref_p;
model.p_global_val = x_target;

% running + terminal cost
Q = diag([0.0,0.0, 1e-1,1e-1]);
Q_f = diag([1e-5,1e-5, 1e1,1e1]);
model.f_q   = 1e-4*norm_2(model.u)^2 + (model.x - model.p_global)'*Q*(model.x - model.p_global);  % small control effort
model.f_q_T = (model.x - model.p_global)'*Q_f*(model.x - model.p_global);


%% Time discretization
problem_options.T = T;
problem_options.N_stages = N_stages;
problem_options.N_finite_elements = N_fe;
problem_options.n_s = 1;
problem_options.cross_comp_mode= "FE_FE";

%% Homotopy / NLP options
reg_solver_options.complementarity_tol = 1e-6;
reg_solver_options.print_level = 0;
reg_solver_options.opts_casadi_nlp.ipopt.max_iter = 500;
% reg_solver_options.homotopy_update_rule = 'superlinear';
% solver_options.homotopy_update_exponent = 2;
reg_solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; % use 'mumps' if needed
%% Solve
ocp = nosnoc.ocp.Solver(model, problem_options, reg_solver_options);
ocp.solve();

%% Create sim model and integrator
if model_plant_mismatch
    sim_model = nosnoc.model.PDSObjects();
    sim_model.u  = u1;
    sim_model.u0 = [0;0];
    sim_model.lbu = [inf;-inf];
    sim_model.ubu = [ inf; inf];

    manip_sim = nosnoc.objects.Ball(R, n_d);      % actuated
    obj_sim   = nosnoc.objects.Ball(R_obj, n_d);  % passive

    % Initial states
    manip_sim.x0 = [-5; 5];
    obj_sim.x0   = [ 0; 0];

    % Dynamics
    manip_sim.f_rhs = u1;        % actuated
    obj_sim.f_rhs   = [0; 0];    % passive

    % (Optional) simple bounds to keep things in view
    manip_sim.lbx = [-inf; -inf];  manip_sim.ubx = [inf; inf];
    obj_sim.lbx   = [-inf; -inf];  obj_sim.ubx   = [inf; inf];

    sim_model.addContact(manip_sim, obj_sim);

    sim_problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
    integrator_options = nosnoc.integrator.Options();
    sim_solver_options = integrator_options.fesd_solver_opts; % Initialize all options related to the MPEC solver used for solving nosonc problems.
    % Choosing the Runge - Kutta Method and number of stages
    sim_problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
    sim_problem_options.n_s = 2; % Number of stage points in the RK method (determines accuracy)

    % Time-settings
    sim_problem_options.N_finite_elements = 4; % Number of finite element (integration steps) on every control interval (optionally a vector might be passed).
    sim_problem_options.T_sim = problem_options.h;
    sim_problem_options.N_sim = 5;
    sim_problem_options.print_level = 0;
    sim_problem_options.cross_comp_mode = "FE_FE";

    % Simulation solver options
    sim_solver_options.print_level = 3;
    sim_solver_options.homotopy_steering_strategy = HomotopySteeringStrategy.DIRECT; % Use the $\ell_{\infty}$ steering strategy
    sim_solver_options.complementarity_tol = 1e-8; % Value to drive the complementarity residual to.
    % sim_solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; % Using an HSL solver is a significant speedup over the default 'mumps' but requires installation
    integrator = nosnoc.Integrator(sim_model, sim_problem_options, integrator_options);
end

%% Extract single MPC
x_res = ocp.get("x");   % states (manip then object)
u_res = ocp.get("u");   % controls
h_res = ocp.get("h");   % step sizes
t_res = ocp.get_time_grid();

%% MPC options
mpc_options.fast_sigma_0 = 1; %
mpc_options.solve_advanced_problem = true;
mpc_options.advanced_n_qpecs = 1;
mpc_options.discard_constraints_in_hessian = true;
mpc_options.do_shift_initialization = false;
mpc_options.warmstart_qpec = false;
mpc_options.advanced_problem_type = 'full';
% problem_options.lift_complementarities = true;

%% create mpc object
if use_rtmpc
    mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, ccopt_options_2, ccopt_options_2);
    % mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, solver_options, mpecopt_options);
    % mpc = nosnoc.mpc.Rtmpc(model, mpc_options, problem_options, mpecopt_options, lcqpow_options);
else
    % mpc = nosnoc.mpc.FullMpcc(model, mpc_options, problem_options, ccopt_options);
    mpc = nosnoc.mpc.FullMpcc(model, mpc_options, problem_options, reg_solver_options);
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
qp_solved = [];
qp_accepted = [];
qpec_solved = [];


for step=1:N_steps
    [u_i, stats] = mpc.get_feedback(x0);
    if use_rtmpc
        qp_solved = [qp_solved, stats.qp_solved];
        qp_accepted = [qp_accepted, stats.qp_accepted];
        qpec_solved = [qpec_solved, stats.qpec_solved];
        if ~stats.qpec_solved
            pause
        end
    else
        qp_solved = [qp_solved, nan];
        qp_accepted = [qp_accepted, nan];
        qpec_solved = [qpec_solved, nan];
    end
    % f_opt = [f_opt, mpc.get_objective()];
    tf_i = stats.feedback_time;
    tf = [tf, tf_i];
    fprintf("MPC step: %d, Feedback time: %d\n", step, tf_i);
    fprintf("Control input: %2.2f\n", u_i);
    if model_plant_mismatch
        integrator.set_x0(x0);
        [t_grid, x_sim] = integrator.simulate("u", repmat(u_i, [1,5]), "x0", x0);
        x0 = x_sim(:, end);
    else
        x0 = mpc.get_predicted_state();
        % if mod(N_steps,6) == 0
        %  x0(1:2) = x0(1:2) + 0.1-0.2*rand(2,1);
        % end
    end
    mpc.do_preparation(x0);

    % Change references!
    if step == 40
        mpc.set_param('p_global',[],[manip_target; [4; 4]]);
    end

    if step == 70
        mpc.set_param('p_global',[],[manip_target; [2; 0]]);
    end

    if step == 90
        mpc.set_param('p_global',[],[manip_target; [0; 0]]);
    end

    x = [x, x0];
    u = [u,u_i];
    t = [t, t(end) + problem_options.h];
    %pause
end

%% Plot all

fig = figure('Position',[50 50 1200 700]); hold on; axis equal
xlim([-10,10]); ylim([-2,10]); box on
set(gca,'XTick',[-10,0,10], 'YTick',[0,5,10], 'LineWidth', 2.5);
xlabel('$c_x$','Interpreter','latex'); ylabel('$c_y$','Interpreter','latex');

% animate with your helper
plot_pass_discs(diff(t), x, [R, R_obj], ["circle","circle"], fig, 'mpc_pds');
title('One actuated disc pushing one passive disc')
