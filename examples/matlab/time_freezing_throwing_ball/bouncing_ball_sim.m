function [tout,yout,error] = bouncing_ball_sim(u_opt,T_opt,N_stages,y0,beta,e,q_target )

tol = 1e-12;
dT = T_opt/N_stages;
tstart = 0;
tfinal = dT;

y_final = q_target;

refine = 4;
options1 = odeset('Events',@events1,'OutputSel',1,'RelTol',tol,'AbsTol',tol/10,'Refine',refine);


tout = tstart;
yout = y0.';
teout = [];
yeout = [];
ieout = [];

mode = 0;
mode_opt = [mode];
bounces = 0;
for i = 1:N_stages
    u = u_opt(:,i);
    % Solve until the first terminal event.
    t = tstart;
    tout_i = [];
    yout_i = [];
    while abs(t(end)-tfinal)>tol
        if abs(tstart-tfinal)>tol
            [t,y,te,ye,ie] = ode15s(@(t,y) ball_dynamics(y,u),[tstart tfinal],y0,options1);
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
                y0(4) =  -e*y0(4);
                bounces = bounces+1;
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
    % Set the new initial conditions
%     y0 = y(end,:);
    tstart = t(end);
    tfinal = tstart + dT;
end
figure
plot(yout(:,1),yout(:,2),'LineWidth',1.5)
hold on
yline(0,'k-','LineWidth',1.5)
plot(y_final(1),y_final(2),'r.','MarkerSize',10)
xlabel('$q_x$','Interpreter','latex')
ylabel('$q_y$','Interpreter','latex')
ylim([min(yout(:,2))-1 max(yout(:,2))+1])
xlim([0 yout(end,1)])
%     grid on
try
    plot(yeout(:,1),yeout(:,2),'rx')
catch
end
error = norm(y_final-yout(end,1:2)');
fprintf('Terminal error is : %2.2e \n',error);
fprintf('Number of bounces: %d \n',bounces);
% --------------------------------------------------------------------------

    function [value,isterminal,direction] = events1(t,y)
        % Locate the time when height passes through zero in a decreasing direction
        % and stop integration.
        value = y(2);     % detect height = 0
        isterminal = 1;   % stop the integration
        direction = -1;   % negative direction
    end
% --------------------------------------------------------------------------
    function [f_y] = ball_dynamics(y,u)
        % Locate the time when height passes through zero in a decreasing direction
        % and stop integration.
        q = y(1:2,:);
        v = y(3:4,:);
        drag = sqrt(v(1)^2^2+v(2)^2+1e-3);
        f_y = [v;[0;-9.81]+u-beta*v*drag];
    end
% --------------------------------------------------------------------------

end


