function [tout,yout,n_bounces] = two_springs_matlab(T_sim,y0,e,tol)

tstart = 0;
tfinal = T_sim;


refine = 4;
options1 = odeset('Events',@events1,'OutputSel',1,'RelTol',tol,'AbsTol',tol/10,'Refine',refine);

tout = tstart;
yout = y0.';
teout = [];
yeout = [];
ieout = [];

mode = 0;
mode_opt = [mode];
n_bounces = 0;

% Solve until the first terminal event.
t = tstart;
tout_i = [];
yout_i = [];
while abs(t(end)-tfinal)>tol
    if abs(tstart-tfinal)>tol
        [t,y,te,ye,ie] = ode45(@(t,y) twospring_dynamics(y),[tstart tfinal],y0,options1);
        teout = [teout; te];          % Events at tstart are never reported.
        yeout = [yeout; ye];
        tout_i = [tout_i ;t];
        yout_i = [yout_i ;y];
    end

    if isequal(ie,1)
        if abs(te-tstart)<tol
            tstart = tfinal;
        else
            mode = 1-mode;
            tstart = te;
            y0 = y(end,:);
            y0(3) =  -e*y0(3);
            n_bounces = n_bounces+1;
        end
        ie = [];
        ye = nan;
    else
        y0 = y(end,:);
    end
end

mode_opt  = [mode_opt;mode];
tout = [tout; tout_i];
yout = [yout; yout_i];
ieout = [ieout; ie];
tstart = t(end);
tfinal = tstart + T_sim;

figure
subplot(311)
plot(tout,yout(:,1),'LineWidth',1.5)
hold on
plot(tout,yout(:,2),'LineWidth',1.5)
grid on
yline(0.2,'k-','LineWidth',1.5)
ylabel('$q$','interpreter','latex');
xlabel('$t$','interpreter','latex');
subplot(312)
plot(tout,yout(:,3),'LineWidth',1.5)
hold on
plot(tout,yout(:,4),'LineWidth',1.5)
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
fprintf('Number of bounces: %d \n',n_bounces);
subplot(313)
stem(tout,[nan;diff(yout(:,3))])
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\Lambda_{\mathrm{n}}$','interpreter','latex');
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


