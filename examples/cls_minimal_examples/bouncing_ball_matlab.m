function [t_grid,x_traj,n_bounces, lambda_normal] = bouncing_ball_matlab(T_sim,x0,e,tol)

tstart = 0;
tfinal = T_sim;


refine = 4;
options1 = odeset('Events',@events1,'OutputSel',1,'RelTol',tol,'AbsTol',tol/10,'Refine',refine);

t_grid = tstart;
x_traj = x0.';
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
        [t,y,te,ye,ie] = ode15s(@(t,y) bouncing_ball_dynamics(y),[tstart tfinal],x0,options1);
        teout = [teout; te];          % Events at tstart are never reported.
        yeout = [yeout; ye];
        t_grid = [t_grid ;t];
        x_traj = [x_traj ;y];
    end

    if isequal(ie,1)
        if abs(te-tstart)<tol
            tstart = tfinal;
        else
            mode = 1-mode;
            tstart = te;
            x0 = y(end,:);
            x0(2) =  -e*x0(2);
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
grid on
ylabel('$q$ [m]','interpreter','latex');
xlabel('$t$ [s]','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

subplot(n_subplot, 1, 2)
plot(t_grid,x_traj(:,2),'LineWidth',1.5)
grid on
xlabel('$t$ [s]','interpreter','latex');
ylabel('$v$ [m/s]','interpreter','latex');
fprintf('Number of bounces: %d \n',n_bounces);
set(gca,'TickLabelInterpreter','latex');


lambda_normal = [0;diff(x_traj(:,2))];
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
        value = y(1);     % detect height = 0
        isterminal = 1;   % stop the integration
        direction = -1;   % negative direction
    end
% --------------------------------------------------------------------------
    function [f_y] = bouncing_ball_dynamics(y,u)
        % Locate the time when height passes through zero in a decreasing direction
        % and stop integration.
        q = y(1,:);
        v = y(2,:);
        m = 1;
        g = 9.81;
        f_v = [-m*g]./m;
        f_y = [v;f_v];
    end
% --------------------------------------------------------------------------

end


function [outputArg1,outputArg2] = untitled(inputArg1,inputArg2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end