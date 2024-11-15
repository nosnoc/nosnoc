%%
clear; clc; close all
import casadi.*
import nosnoc.*
%%
plot_integrator_output = 1;
plot_continious_time_sol = 1;
%% discretization settings
T_sim = pi/2;
N_sim  = 29;
N_finite_elements = 2;
R_osc  = 1;

%% Options and model objects
problem_options = nosnoc.Options();
solver_options = nosnoc.solver.Options();
model = nosnoc.model.Pss();
%% settings
problem_options.rk_scheme = RKSchemes.RADAU_IIA; %'Gauss-Legendre';
problem_options.n_s = 4;
problem_options.dcs_mode = 'Heaviside'; % 'Step;
problem_options.N_finite_elements = N_finite_elements;
problem_options.T_sim = T_sim;
problem_options.N_sim = N_sim;
problem_options.print_level  = 1;

%% Time settings
x_star = [exp(1);0];
T = T_sim;
x_star = [exp(T-1)*cos(2*pi*(T-1));-exp((T-1))*sin(2*pi*(T-1))];

omega = -2*pi;
A1 = [1 omega;...
    -omega 1];
A2 = [1 -omega;...
    omega 1];
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1;x2];
c = x'*x-1;
f_1 = A1*x;
f_2 = A2*x;
F = [f_1 f_2];
x0 = [exp(-1);0];
model.x0 = x0;

model.x = x;
model.c = c;
model.S = [-1;1];
model.F = F;
%% Call integrator
integrator = nosnoc.integrator.FESD(model, problem_options, solver_options);
[t_grid, x_res] = integrator.simulate();
%% numerical error
x_fesd = x_res(:,end);
error_x = norm(x_fesd-x_star,"inf");
fprintf(['Numerical error with h = %2.3f and ' char(problem_options.rk_scheme) ' with n_s = %d stages is: %5.2e: \n'],problem_options.h_sim,problem_options.n_s,error_x);
%% plot_solution_trajectory
t_star = R_osc; % eact switching time
h_opt_full = integrator.get('h');
x1_opt = x_res(1,:);
x2_opt = x_res(2,:);
%%
if plot_integrator_output
    figure
    subplot(121)
    plot(t_grid,x1_opt,'linewidt',1.0);
    grid on
    hold on
    plot(t_grid,x2_opt,'linewidt',1.0);
    hh = -3:1:3;
    plot(hh*0+t_star,hh,'k')
    xlabel('$t$','interpreter','latex');
    ylabel('$x(t)$','interpreter','latex');
    legend({'$x_1(t)$','$x_2(t)$'},'interpreter','latex');
    subplot(122)
    plot(x1_opt,x2_opt,'linewidt',1.8);
    hold on
    theta = 0:0.01:2*pi;
    x = R_osc*(cos(theta));
    y = R_osc*(sin(theta));
    plot(x,y,'r','linewidth',1.5)
    xlabel('$x_1$','interpreter','latex');
    ylabel('$x_2$','interpreter','latex');
    axis equal
    grid on
    [X,Y] = meshgrid(-3:0.35:3,-3:0.35:3);
    [U,V] = oscilator_eval(X,Y);
    quiver(X,Y,U,V,'Color',0.65*ones(3,1),'LineWidth',1.2,'AutoScaleFactor',3);
    xlim([-exp(1) exp(1)]*1.01)
    ylim([-exp(1) exp(1)]*0.9)
    %
    figure
    stairs(h_opt_full)
    ylim([min(h_opt_full)*0.8 max(h_opt_full)*1.2])
    xlabel('integration step')
    ylabel('$h$','Interpreter','latex');
end

%% plot_continious_time_sol
x_res_extended = integrator.get_full('x');
tgrid_long  = integrator.get_time_grid_full();
%%
if plot_continious_time_sol
    % tgrid_long = results.extended.t_grid;

    %
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
        p1 = polyfit(tt,xx1,length(xx1)-2);
        p2 = polyfit(tt,xx2,length(xx2)-2);
        t_eval = linspace(tt(1),tt(end),50);
        yy1 = polyval(p1,t_eval);
        yy2 = polyval(p2,t_eval);
        % store
        x1_very_fine = [x1_very_fine,yy1];
        x2_very_fine = [x2_very_fine,yy2];
        tgrid_very_fine = [tgrid_very_fine,t_eval];
        % plot
        plot(t_eval,yy1,'b');
        hold on;
        plot(t_eval,yy2,'r');
        plot(tt,xx1,'b.');
        plot(tt,xx2,'r.');
        grid on
    end
    xline(R_osc,'m')
    for ii = 1:length(t_grid)
        xline(t_grid(ii),'k--')
    end
end
