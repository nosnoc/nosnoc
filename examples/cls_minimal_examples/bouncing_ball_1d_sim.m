clear all;
% clc;
import casadi.*
close all
%% init nosnoc settings and model
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
model = nosnoc.model.Cls();
%% Simulation setings
N_FE = 2;
T_sim = 3;
N_sim = 20;

problem_options.T_sim = T_sim;
problem_options.N_sim = N_sim;
problem_options.N_finite_elements = N_FE;
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
problem_options.n_s = 2;
problem_options.cross_comp_mode = 1;
problem_options.dcs_mode = DcsMode.CLS;
problem_options.no_initial_impacts = 1;

% Initialization 
problem_options.initial_Lambda_normal = 0;
problem_options.initial_lambda_normal = 0;
problem_options.initial_Y_gap = 0;
problem_options.initial_y_gap = 0;

%problem_options.fixed_eps_cls = 1;
%problem_options.relax_terminal_numerical_time = ConstraintRelaxationMode.ELL_1;
%problem_options.rho_terminal_numerical_time = 1e3;
%problem_options.relax_fesdj_impulse = ConstraintRelaxationMode.ELL_2;
%problem_options.rho_fesdj_impulse = 1e6;
%problem_options.gamma_h = 0.99;

solver_options.decreasing_s_elastic_upper_bound = 1; % elasic mode with decreasing bounds for the elstaic slacks
solver_options.print_level = 2;


%% model defintion
g = 9.81;
x0 = [0.8;0];

q = SX.sym('q',1);
v = SX.sym('v',1);
model.M = 1;
model.x = [q;v];
model.e = 1;
model.mu = 0;
model.x0 = x0;
model.f_v = -g;
model.f_c = q;

%% Analytic solution

if model.e == 0
    t_s = sqrt(2*x0(1)/g);
    tt1 = linspace(0,t_s,100);
    tt2 = linspace(t_s,T_sim,100);
    v1 = x0(2)-g*tt1;
    q1 = x0(1)+x0(2)*tt1-g*tt1.^2/2;
    if model.e == 0
        v2 = 0*(tt2-t_s);
        q2 = 0*(tt2-t_s);
    else
        v2 = -model.e*v1(end)-g*(tt2-t_s);
        q2 = q1(end)+v2(1)*(tt2-t_s)-g*(tt2-t_s).^2/2;
    end
    a_t_grid = [tt1,tt2];
    q_traj= [q1,q2];
    v_traj = [v1,v2];
    x_traj = [q_traj',v_traj'];
    lambda_normal = [q1*0,q2*0+g];
    Lambda_star = v2(1)-v1(end);
else
    [a_t_grid,x_traj,n_bounces, lambda_normal] = bouncing_ball_matlab(T_sim,model.x0,model.e,1e-12);
    Lambda_star = max(abs(diff(x_traj(:,2))));

end

%% Call nosnoc integrator
integrator = nosnoc.Integrator(model, problem_options, solver_options);
[t_grid, x_res, t_grid_full, x_res_full] = integrator.simulate();

%% read and plot results
qx = x_res(1,:);
vx = x_res(2,:);
Lambda_normal = integrator.get("Lambda_normal");

figure
subplot(311)
plot(t_grid,qx);
hold on
plot(a_t_grid,x_traj(:,1))
legend({'$q$ - numerical','$q$ - analytic'},'interpreter','latex');
xlim([0 t_grid(end)])
ylim([-0.1 x0(1)+0.1])
grid on
ylabel('$q$','interpreter','latex');
xlabel('$t$','interpreter','latex');
subplot(312)
plot(t_grid,vx);
hold on
plot(a_t_grid,x_traj(:,2))
legend({'$q$ - numerical','$q$ - analytic'},'interpreter','latex');
ylim([-3 3])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
subplot(313)
if problem_options.no_initial_impacts == 1
    Lm = [];
    for ii=0:length(Lambda_normal)-1
        if mod(ii,N_FE-1) == 0
            Lm = [Lm,nan,nan];
        end
        Lm = [Lm, Lambda_normal(ii+1)];
    end
    stem(t_grid,[Lm,nan])
else
    stem(t_grid,[Lambda_normal,nan])
end
hold on
yline(Lambda_star,'k--')
xlim([-0.01 t_grid(end)])
ylim([-0.1 max([Lambda_normal,Lambda_star])+1])
grid on
legend({'$\Lambda$ - numerical','$\Lambda$ - analytic'},'interpreter','latex');
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda$','interpreter','latex');

fprintf('position error %2.2e \n',abs(x_traj(end,1)-qx(end)))
fprintf('velocity error %2.2e \n',abs(x_traj(end,2)-vx(end)))

