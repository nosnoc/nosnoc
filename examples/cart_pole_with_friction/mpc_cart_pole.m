clear; clc; close all;
import casadi.*
%% plot and video settings
latexify_plot()
nice_plot_colors;
%%
filename_pdf = 'cartpole_mpc';
video_name = filename_pdf;
linewidth = 1.5;
save_video = 0;
create_video = 0;
real_time_mpc_plot = 1;
use_mpecopt = 0;
%% Model settings
% Discretization options
model_plant_missmatch = 0;
model_plant_same_integrator = 1;
N_sim = 1;

N_steps = 40; % number of MPC steps in simulation;
N_stages = 15; % number of control stages;

DT = 0.2;
F_friction_model = 2; % Friction force amplitude
F_friction_plant = 2;

N_FE = 2;
use_fesd = 1;
gamma_h = 1;

%% solver options
tol = 1e-6;
N_homotopy = 6;
sigma_0 = 1;
fullmpcc_fast_sigma_0 = 1e-3; % sigma for mpec warmstart in mpc

%% nosnoc setup
problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
mpc_options = nosnoc.mpc.Options(); % Initialize all options related to the MPC implementation.
solver_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosonc problems.
%% Optimization solver 

mpecopt_options = mpecopt.Options();
mpecopt_options.compute_tnlp_stationary_point = false;
mpecopt_options.settings_lpec.lpec_solver = "Gurobi";
mpecopt_options.settings_casadi_nlp.ipopt.linear_solver = 'ma27';
mpecopt_options.rho_TR_phase_i_init = 1e1;
mpecopt_options.rho_TR_phase_ii_init = 1e-2;

%% Model variables and parameters
% specify initial and desired state
x0 = [1; 0/180*pi; 0; 0]; % start downwards
x_ref = [0; 180/180*pi; 0; 0]; % end upwards
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
q = vertcat(px, theta);
% Velocities
v = SX.sym('v');
theta_dot = SX.sym('theta_dot');
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
Cor = [0, -m2 * link_length*theta_dot*sin(theta);...
    0,   0];
f_all = [0; -m2*g*link_length*sin(theta)] + [u; 0] - Cor*q_dot;

% Dynamics for v > 0
f_1 = [q_dot;...
    inv(M)*(f_all-[F_friction; 0])];
% Dynamics for v<0
f_2 = [q_dot;...
    inv(M)*(f_all+[F_friction; 0])];

f_1_fun = Function('f_1_fun',{x,u,F_friction},{f_1});
f_2_fun = Function('f_2_fun',{x,u,F_friction},{f_2});

F = [f_1_fun(x,u,F_friction_model), f_2_fun(x,u,F_friction_model)];

ubx = [5; inf; inf; inf];
lbx = [-5; -inf; -inf; -inf];
u_max = 30;

Q = diag([100; 10; 1; 1]);
R = 1;
% stage cost
f_q = (x-y_ref_p(1:4))'*Q*(x-y_ref_p(1:4))+ (u-y_ref_p(5))'*R*(u-y_ref_p(5)); % running/stage costs
% terminal cost
Q_terminal = Q*25;
f_q_T  = (x-y_ref_p(1:4))'*Q_terminal*(x-y_ref_p(1:4));
% Populate nosnoc model
model.c = v;         % switching function c: cart velocity
model.S = [1; -1];   % sign matrix S % f_1 for c>0, f_2 for c<0
model.p_global = y_ref_p;
model.p_global_val = [x_ref;0];
model.F = F;
model.f_q = f_q;
model.f_q_T = f_q_T;
model.lbx = lbx;
model.ubx = ubx;
model.x = x;
model.x0 = x0;
model.u = u;
model.lbu = -u_max;
model.ubu = u_max;

% model.g_terminal = model.x-x_ref;
% problem_options.relax_terminal_constraint = "ELL_INF";

%% OCP problem options
problem_options.T = N_stages*DT;  % Time horizon
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 2;
problem_options.dcs_mode = 'Stewart';
problem_options.N_stages = N_stages; % number of control intervals
problem_options.N_finite_elements = N_FE; % number of finite element on every control interval
problem_options.cross_comp_mode = "FE_FE";
problem_options.use_fesd = use_fesd;
problem_options.gamma_h = gamma_h;


%% MPEC solver options 
solver_options.complementarity_tol = tol; % Value to drive the complementarity residual to.
solver_options.N_homotopy = N_homotopy; % Maximum number of homotopy iterations.
solver_options.print_level = 3;
solver_options.sigma_0 = sigma_0;
solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; %
solver_options.homotopy_steering_strategy = "ELL_INF";

%% MPC options
mpc_options.fullmpcc_fast_sigma_0 = fullmpcc_fast_sigma_0; % Set the fast sigma for followup solves to the complementarity tol to only do 1 nlp solve.

%% create plant
if model_plant_missmatch
    % Integrator model
    sim_model = nosnoc.model.Pss();
    sim_model.c = v;  
    sim_model.x = x;
    sim_model.x0 =  x0;
    sim_model.u = u;
    sim_model = get_cart_pole_with_friction_model(true, F_friction_plant);
    sim_model.lbx = -inf(n_x,1); sim_model.ubx = inf(n_x,1); 
    sim_model.lbu = -inf;  sim_model.ubu = inf;
    % Dynamics of the piecewise smooth systems
    sim_problem_options = nosnoc.Options(); % Initialize all options related to the optimal control or simulation problem.
    % sim_solver_options = nosnoc.reg_homotopy.Options(); % Initialize all options related to the MPEC solver used for solving nosonc problems.
    integrator_options = nosnoc.integrator.Options();
    integrator_options.integrator_plugin = "FESD";
    

    sim_solver_options = integrator_options.fesd_solver_opts; % the fesd integrator uses an mpec solver, call and modify its options

    % Choosing the Runge - Kutta Method and number of stages
    sim_problem_options.rk_scheme = RKSchemes.RADAU_IIA; % Type of scheme
    if model_plant_same_integrator
        sim_problem_options.n_s = 2; % Number of stage points in the RK method (determines accuracy)
        sim_problem_options.N_finite_elements = N_FE; % Number of finite element (integration steps) on every control interval (optionally a vector might be passed).
        sim_problem_options.N_sim = 1;
        N_sim = 1;
    else
        sim_problem_options.n_s = 3; % Number of stage points in the RK method (determines accuracy)
        sim_problem_options.N_finite_elements = 2; % Number of finite element (integration steps) on every control interval (optionally a vector might be passed).
        sim_problem_options.N_sim = N_sim;
    end
    sim_problem_options.T_sim = DT;
    sim_solver_options.print_level = 0;
    sim_solver_options.homotopy_steering_strategy = HomotopySteeringStrategy.ELL_INF; % Use the $\ell_{\infty}$ steering strategy
    sim_solver_options.complementarity_tol = 1e-7; % Value to drive the complementarity residual to.
    sim_solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27'; % Using an HSL solver is a significant speedup over the default 'mumps' but requires installation

    % create nosnoc integrator
    integrator = nosnoc.Integrator(sim_model, sim_problem_options, integrator_options);
end

if real_time_mpc_plot && model_plant_missmatch
    figure('Position', [100, 100, 2000, 1000]); % [left, bottom, width, height]
end
if save_video && real_time_mpc_plot
    outputVideo = VideoWriter([video_name '.mp4'], 'MPEG-4');
    outputVideo.FrameRate = 5; % Adjust the frame rate as needed
    open(outputVideo);
end

% create mpc object
if use_mpecopt
    mpc_options.fullmpcc_do_shift_initialization = false;
    mpc = nosnoc.mpc.BnlpLpcc(model, mpc_options, problem_options, mpecopt_options);
else
    mpc_options.fullmpcc_do_shift_initialization = true;
    mpc = nosnoc.mpc.FullMpcc(model, mpc_options, problem_options, solver_options);
    % mpc = nosnoc.mpc.FullMpcc(model, mpc_options, problem_options, mpecopt_options);
    % mpecopt_options
end



x0 = model.x0;
u = [];
t = 0;
tf = [];
x = x0;
d_hat = zeros(4,1);
% error_in_x0 = 0;
error_in_x0 = zeros(4,1);
error_in_predicted_states = [];
f_opt = [];
x_high_res = x0;
t_high_res = 0;
x_estimator_pred = x0;
x_estimator_true = x0;
x_mpc_pred = x0;
t_estimator_actual = 0;
t_estimator_pred = 0;
e_over_time = zeros(n_x,1);
dynamic_targ = [x_ref;0];
if model_plant_missmatch
    for step=1:N_steps
        % Estimate state and distrubnace
        % y = C*x0; % y = C*x_high_res(:,end); % true state (take measurment)
        % % error between measurment and state estimator prediction
        % e = y - C*x_estimator_pred(:,end) - Cd*d_estimator_pred(:,end);
        % % correct the state estimator;
        % e_over_time = [e_over_time,e];
        % x_estimator_true = [x_estimator_true, x_estimator_pred(:,end)+Lx*e];
        % d_estimator_true = [d_estimator_true, d_estimator_pred(:,end)+Ld*e];
        % t_estimator_actual = [t_estimator_actual, t_estimator_actual(end)+DT];
        [u_i, stats] = mpc.get_feedback(x0);
        % f_opt = [f_opt, mpc.get_objective()];
        tf_i = stats.feedback_time;
        tf = [tf, tf_i];
        fprintf("MPC step: %d, Feedback time: %d\n", step, tf_i);
        fprintf("Control input: %2.2f\n", u_i);
        mpc.do_preparation();

        % what would the plant to? - compute model plat missmatch for current u_mpc
        if real_time_mpc_plot
            x0_local = x0;
            x_mpc_sim = x0;
            t_grid = mpc.get_time_grid();
            t_grid_u = mpc.get_control_grid();
            u_traj = mpc.get('u');
            x_traj = mpc.get('x');
            for jj = 1:N_stages
                ind_fe = (jj-1)*N_FE+1:(jj)*N_FE+1;
                % get results
                
                tspan_jj = t_grid (ind_fe)-t_grid(ind_fe(1));
                integrator.set_x0(x0_local);
                [t_mpc_jj, x_mpc_jj] = integrator.simulate("u", repmat(u_traj(jj), [1,N_sim]), "x0", x0_local);
                if model_plant_same_integrator
                    x_mpc_sim = [x_mpc_sim, x_mpc_jj(:,2:end)];
                else
                    x_mpc_sim = [x_mpc_sim, x_mpc_jj(:,[3,5])];
                end
                x0_local = x_mpc_jj(:,end);
            end
            % error_in_predicted_controls = vecnorm(stats.x-x_mpc);
            error_in_predicted_states = abs(x_traj-x_mpc_sim);
        end
        % Advance plat in time
        integrator.set_x0(x0);
        [t_grid_sim, x_sim] = integrator.simulate("u", repmat(u_i, [1,N_sim]), "x0", x0);
        % Store result here:
        % error_in_x0 = [error_in_x0,norm(x_sim(:,end)-stats.x(:,N_FE+1))];
        % error_in_x0 = [error_in_x0,x_sim(:,end)-stats.x(:,N_FE+1)];
        x_high_res = [x_high_res,x_sim];
        t_high_res = [t_high_res, t_high_res(end) + t_grid_sim];
        change_in_state  = norm(x0-x_sim(:, end));
        x0 = x_sim(:, end);
        x = [x, x0];
        u = [u,u_i];
        t = [t, t(end) + problem_options.h];
        fprintf("Change in state: %2.2e\n",change_in_state);

        % MPC real time plot
        if real_time_mpc_plot
            if ~save_video
                clf
            end
            % position
            subplot(321)           
            plot(t_high_res,x_high_res(1,:),'LineWidth',linewidth,'Color',matlab_blue);
            hold on;
            plot(t_grid+t_high_res(end-length(t_grid_sim)),x_traj(1,:),'LineWidth',linewidth,'Color', matlab_red);
            xlabel("$t$")
            ylabel("$q$")
            grid on
            yline(0,'k--','LineWidth',1.5)
            xline(t_high_res(end-length(t_grid_sim)),'k-')
            xline(t_high_res(end),'k-')
            xlim([0 N_steps*DT]);
            ylim([-2 2])
            % angle
            subplot(323)
            plot(t_high_res,x_high_res(2,:),'LineWidth',linewidth,'Color',matlab_blue);
            hold on
            plot(t_grid+t_high_res(end-length(t_grid_sim)),x_traj(2,:),'LineWidth',linewidth,'Color', matlab_red);
            xlabel("$t$")
            ylabel("$\theta$")
            ylim([-pi*1.1 pi*1.1])
            yline(pi,'k--','LineWidth',linewidth-0.5)
            xline(t_high_res(end-length(t_grid_sim)),'k-')
            xline(t_high_res(end),'k-')
            xlim([0 N_steps*DT]);
            grid on
            subplot(325)
            stairs(t,[u,u(end)],'LineWidth',linewidth,'Color',matlab_blue);
            hold on
            stairs(t_grid_u+t(end-1),[u_traj,u_traj(end)],'LineWidth',linewidth,'LineStyle','-','Color',matlab_red);
            xlabel("$t$")
            ylabel("$u$")
            xlim([0 N_steps*DT]);
            grid on
            u_max = 30;
            ylim([-1.1*u_max 1.1*u_max])
            hold on
            yline(-30,'k--')
            yline(30,'k--')
            yline(-F_friction_plant,'k-')
            yline(F_friction_plant,'k-')
            xline(t(end-1),'k-')
            xline(t(end),'k-')
            % distance to ref in position
            subplot(322)
            semilogy(t_high_res,abs(x_high_res(1,:)-x_ref(1)),'LineWidth',linewidth,'Color',matlab_blue);
            hold on;
            semilogy(t_grid+t_high_res(end-length(t_grid_sim)),abs(x_traj(1,:)-x_ref(1)),'LineWidth',linewidth,'Color', matlab_red);
            xlabel("$t$")
            ylabel("$|q-q_{|\mathrm{ref}}|$")
            xline(t_high_res(end-length(t_grid_sim)),'k-')
            xline(t_high_res(end),'k-')
            grid on
            % yline(1e-8,'k--','LineWidth',1.5)
            ylim([1e-7 1e1])
            xlim([0 N_steps*DT]);

            subplot(324)
            semilogy(t_high_res,abs(x_high_res(2,:)-x_ref(2)),'LineWidth',linewidth,'Color',matlab_blue);
            hold on;
            semilogy(t_grid+t_high_res(end-length(t_grid_sim)),abs(x_traj(2,:)-x_ref(2)),'LineWidth',linewidth,'Color', matlab_red);
            xlabel("$t$")
            ylabel("$|\theta-\theta_{|\mathrm{ref}}|$")
            grid on
            xline(t_high_res(end-length(t_grid_sim)),'k-')
            xline(t_high_res(end),'k-')
            % ylim([-1 1]);
            ylim([-0.1 pi*1.1])
            % yline(1e-8,'k--','LineWidth',linewidth-0.5)
            ylim([1e-7 1e1])
            xlim([0 N_steps*DT]);
            grid on
            hh = subplot(326);
            % semilogy(t,error_in_x0,'LineWidth',linewidth,'Color',matlab_blue);
            % plot(t,error_in_x0,'LineWidth',linewidth);
            % hold on
            % semilogy(t(end-1)+stats.t_grid,error_in_predicted_controls+1e-16,'LineWidth',linewidth,'Color',matlab_red);
            plot(t(end-1)+t_grid,error_in_predicted_states,'LineWidth',linewidth);
            % error_in_predicted_controls(1)
            hold on
            xlabel("$t$")
            % ylabel("$|x^\mathrm{model}(t;u^\mathrm{mpc}) - x^\mathrm{plant}(t;u^\mathrm{mpc})|$");
            ylabel("$|x_{\mathrm{mpc}}-x_{\mathrm{sim}}|$");
            grid on
            ylim([-3 3])
            hold on
            xline(t(end-1),'k-')
            xline(t(end),'k-')
            xlim([0 N_steps*DT]);
            axesHandles = findall(gcf, 'type', 'axes');
            % set(axesHandles, 'FontSize', 14);
            if save_video
                drawnow;
                % Capture the current frame
                frame = getframe(gcf);
                writeVideo(outputVideo, frame);
                clf;
            end
        end
    end
else
    for step=1:N_steps
        [u_i, stats] = mpc.get_feedback(x0);
        tf_i = stats.feedback_time;
        tf = [tf, tf_i];
        fprintf("MPC step: %d, Feedback time: %d\n", step, tf_i);
        fprintf("Control input: %2.2f\n", u_i);
        fprintf("Offset: %2.2e\n", norm(x0-x_ref));
        mpc.do_preparation();
        x0 = mpc.get_predicted_state();
        x = [x, x0];
        u = [u,u_i];
        t = [t, t(end) + problem_options.h];
        x_high_res = x;
        t_high_res = t;
        t_grid = t(end);
    end
end
if save_video && real_time_mpc_plot
    close(outputVideo);
end

%% Plot
f = figure;
subplot(311)
plot(t_high_res,x_high_res(1,:),'LineWidth',linewidth)
xlabel("$t$")
ylabel("$q$")
grid on
yline(0,'k--','LineWidth',1.5)
subplot(312)
plot(t_high_res,x_high_res(2,:),'LineWidth',linewidth )
xlabel("$t$")
ylabel("$\theta$")
ylim([-0.1 pi*1.1])
yline(pi,'k--','LineWidth',linewidth-0.5)
grid on
subplot(313)
stairs(t,[u,u(end)],'LineWidth',linewidth)
xlabel("$t$")
ylabel("$u$")
grid on
u_max = 30;
ylim([-1.1*u_max 1.1*u_max])
hold on
yline(-30,'r--')
yline(30,'r--')
yline(-F_friction_plant,'r-')
yline(F_friction_plant,'r-')
exportgraphics(f, [ filename_pdf '1.pdf']);

%%
figure
subplot(121)
semilogy(t,vecnorm(x-x_ref));
xlabel("$t$")
ylabel("$|x-x_{\mathrm{sp}}|$")
grid on
subplot(122)
semilogy(t,[nan,f_opt])
xlabel("$t$")
ylabel("$f^*$ (objective of every mpc subproblem)")
grid on
%% Animation
if create_video
    outputVideo = VideoWriter([video_name '_cart.mp4'], 'MPEG-4');
    outputVideo.FrameRate = 10; % Adjust the frame rate as needed
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
        plot(pole_X(1:ii,2),pole_Y(1:ii,2),'color',[1 0 0 0.5],'LineWidth',0.5);
        % cart
        xp = [q1_opt(ii)-cart_width1/2 q1_opt(ii)+cart_height/2 q1_opt(ii)+cart_height/2 q1_opt(ii)-cart_width1/2];
        yp = [0 0 cart_height  cart_height];
        patch(xp,yp,'k','FaceAlpha',0.8)

        % targent
        % pole
        plot([x_ref(1),x_ref(1)+(link_length)*cos(x_ref(2)-pi/2)],...
            [cart_center+0,cart_center+link_length*sin(x_ref(2)-pi/2)],'color',[0 0 0 0.1],'LineWidth',3);
        % cart
        xp = [x_ref(1)-cart_width1/2 x_ref(1)+cart_height/2 x_ref(1)+cart_height/2 x_ref(1)-cart_width1/2];
        yp = [0 0 cart_height  cart_height];
        patch(xp,yp,'k','FaceAlpha',0.1)
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

        % frame = getframe(1);
        % im = frame2im(frame);
        % [imind,cm] = rgb2ind(im,256);
        % if ii == 1
        %     imwrite(imind,cm,filename,'gif', 'Loopcount', inf,'DelayTime', time_step);
        % else
        %     imwrite(imind,cm,filename,'gif', 'WriteMode', 'append','DelayTime', time_step);
        % end
        %
        % if ii<length(q1_opt)
        %     clf;
        % end
    end
    close(outputVideo);
    % exportgraphics(gca, [filename_pdf '2.pdf']);
end


