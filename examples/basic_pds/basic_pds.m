clear all
clc
close all
import casadi.*
plot_continious_time_sol = 1;
plot_integrator = 0;
%% discretization parameters
N_sim = 1;
T_sim = 11*pi/12+sqrt(3);

%% NOSNOC settings
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
problem_options.n_s = 4;
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
%problem_options.rk_representation= RKRepresentation.differential_lift_x; 
problem_options.rk_representation = RKRepresentation.integral;
problem_options.cross_comp_mode = CrossCompMode.STAGE_STAGE;
problem_options.use_fesd = false;
%problem_options.cross_comp_mode = CrossCompMode.FE_FE;
problem_options.N_sim = N_sim;
problem_options.N_finite_elements = 3;
problem_options.T_sim = T_sim;
problem_options.rho_h = 1e-4;
problem_options.gcs_lift_gap_functions = true;
problem_options.step_equilibration = StepEquilibrationMode.l2_relaxed_scaled;
%solver_options.homotopy_steering_strategy = 'ELL_INF';
solver_options.complementarity_tol = 1e-10;
solver_options.print_level = 3;

model = nosnoc.model.Pds();

model.x0 = [sqrt(2)/2;sqrt(2)/2];
x = SX.sym('x',2);
model.x = x;
model.c = x(2)+1;
model.f_x_unconstrained = [x(2);-x(1)];

model.x0 = [sqrt(2);sqrt(2)];

integrator = nosnoc.Integrator(model, problem_options, solver_options);
[t_grid, x_res] = integrator.simulate();
%
if plot_integrator
    figure
    plot(x_res(1,:), x_res(2,:))
    grid on
    xlabel('$x_1$','Interpreter','latex')
    ylabel('$x_2$','Interpreter','latex')
    grid on

    c_fun = casadi.Function('c', {model.x}, {model.c});
    c = full(c_fun(integrator.get_full('x')))';
    x_full = integrator.get_full('x');
    lambda = integrator.get('lambda');

    % extened values over all rk stage points
    t_grid_full = integrator.get_time_grid_full();
    lambda_full = integrator.get_full('lambda');

    %

    figure
    plot(t_grid,lambda,'LineWidth',2)
    hold on
    plot(t_grid,x_res,'LineWidth',1.5)
    grid on
    xlabel('$t$','Interpreter','latex')
    legend({'$\lambda(t)$','$x_1(t)$','$x_2(t)$'},'interpreter','latex');
end
%%
x_res_extended = integrator.get_full('x');
tgrid_long  = integrator.get_time_grid_full();
%%
if plot_continious_time_sol
    lw = 5;
    x1_very_fine = [];
    x2_very_fine = [];
    tgrid_very_fine = [];
    figure
    for ii =  1:problem_options.N_stages*problem_options.N_finite_elements*N_sim
        % read
        ind_now = 1+(ii-1)*(problem_options.n_s):(ii)*(problem_options.n_s)+1;
        tt = tgrid_long(ind_now);
        xx1 = x_res_extended(1,ind_now);
        xx2 = x_res_extended(2,ind_now);
        % fit
        p1 = polyfit(tt,xx1,length(xx1)-1);
        p2 = polyfit(tt,xx2,length(xx2)-1);
        t_eval = linspace(tt(1),tt(end),50);
        yy1 = polyval(p1,t_eval);
        yy2 = polyval(p2,t_eval);
        % store
        x1_very_fine = [x1_very_fine,yy1];
        x2_very_fine = [x2_very_fine,yy2];
        tgrid_very_fine = [tgrid_very_fine,t_eval];
        % plot
        plot(t_eval,yy1,'b', 'Linewidth', lw-3);
        hold on;
        plot(t_eval,yy2,'r', 'Linewidth', lw-3);
        plot(tt,xx1,'b.', 'Markersize', 40);
        plot(tt,xx2,'r.', 'Markersize', 40);
        grid on
    end
    %xline(R_osc,'m')
    for ii = 1:length(t_grid)
        xline(t_grid(ii),'k--')
    end

    % plot closed form
    t1 = 0:0.001:(5*pi/12);
    t2 = (5*pi/12):0.001:(5*pi/12 + sqrt(3));
    t3 = (5*pi/12 + sqrt(3)):0.001:T_sim;

    hold on
    plot(t1,2*sin(t1+pi/4),'Color', [0 0 1 0.4], 'Linewidth', lw);
    plot(t1,2*cos(t1+pi/4),'Color', [1 0 0 0.4], 'Linewidth', lw);
    plot(t2,2*sin(5*pi/12+pi/4) - (t2-5*pi/12),'Color', [0 0 1 0.4], 'Linewidth', lw);
    plot(t2,-1*ones(size(t2)),'Color', [1 0 0 0.4], 'Linewidth', lw);
    plot(t3,sin(t3-(5*pi/12 + sqrt(3))+pi),'Color', [0 0 1 0.4], 'Linewidth', lw);
    plot(t3,cos(t3-(5*pi/12 + sqrt(3))+pi),'Color', [1 0 0 0.4], 'Linewidth', lw);
    hold off
end
