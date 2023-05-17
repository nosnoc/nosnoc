function [t_grid,x_traj,n_bounces, lambda_normal] = two_balls_spring_matlab(T_sim,x0,e,tol)

tstart = 0;
tfinal = T_sim;


refine = 4;
options1 = odeset('Events',@events1,'OutputSel',1,'RelTol',tol,'AbsTol',tol/10,'Refine',refine);

t_grid = [];
x_traj = [];
teout = [];
yeout = [];
ieout = [];

mode = 0;
mode_opt = mode;
n_bounces = 0;

% Solve until the first terminal event.
t = tstart;
while abs(t(end)-tfinal)>tol
    if abs(tstart-tfinal)>tol
        [t,y,te,ye,ie] = ode15s(@(t,y) twospring_dynamics(y),[tstart tfinal],x0,options1);
        teout = [teout; te];          % Events at tstart are never reported.
        yeout = [yeout; ye];
        t_grid = [t_grid ;t];
        x_traj = [x_traj ;y];
    end

    if isequal(ie,1)
        % strip pre impact
        t_grid(end) = [];
        x_traj(end, :) = [];

        if abs(te-tstart)<tol
            tstart = tfinal;
        else
            mode = 1-mode;
            tstart = te;
            x0 = y(end,:);
            x0(3) =  -e*x0(3);
            n_bounces = n_bounces+1;
        end
        ie = [];
        ye = nan;
    else
        x0 = y(end,:);
    end
end

mode_opt  = [mode_opt;mode];
ieout = [ieout; ie];

PLOT_LAMBDA = 0;
n_subplot = 2;
if PLOT_LAMBDA
    n_subplot = n_subplot+1;
end

figure
subplot(n_subplot, 1, 1)
plot(t_grid,x_traj(:,1),'LineWidth',1.5)
hold on
plot(t_grid,x_traj(:,2),'LineWidth',1.5)
grid on
yline(0.2,'k-','LineWidth',1.5)
ylabel('$q$ [m]','interpreter','latex');
xlabel('$t$ [s]','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

subplot(n_subplot, 1, 2)
plot(t_grid,x_traj(:,3),'LineWidth',1.5)
hold on
plot(t_grid,x_traj(:,4),'LineWidth',1.5)
grid on
xlabel('$t$ [s]','interpreter','latex');
ylabel('$v$ [m/s]','interpreter','latex');
fprintf('Number of bounces: %d \n',n_bounces);
set(gca,'TickLabelInterpreter','latex');
legend('ball 1', 'ball 2', 'interpreter','latex')


lambda_normal = [0;diff(x_traj(:,3))];
lambda_normal(abs(lambda_normal) < 2) = 0;
if PLOT_LAMBDA
    subplot(n_subplot, 1, 3)
    stem(t_grid,lambda_normal)
    grid on
    xlabel('$t$','interpreter','latex');
    ylabel('$\Lambda_{\mathrm{n}}$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
end

% --------------------------------------------------------------------------

    function [value,isterminal,direction] = events1(t,y)
        % Locate the time when height passes through zero in a decreasing direction
        % and stop integration.
        R = 0.2;
        value = y(1)-R;     % detect height = 0
        isterminal = 1;   % stop the integration
        direction = -1;   % negative direction
    end
% --------------------------------------------------------------------------
    function [f_y] = twospring_dynamics(y,u)
        % Locate the time when height passes through zero in a decreasing direction
        % and stop integration.
        q = y(1:2,:);
        v = y(3:4,:);
        m = 1;
        l = 1;
        k = 1e4;
        R = 0.2;
        g = 9.81;
        f_v = [-m*g+k*(q(2)-q(1)-l);-m*g-k*(q(2)-q(1)-l)]./m;
        f_y = [v;f_v];
    end
% --------------------------------------------------------------------------

end


