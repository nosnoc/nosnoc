function [tout,yout,error] = car_turbo_sim(u_opt,T_opt,N_stages,hybrid_dynamics)
tol = 1e-10;
dT = T_opt/N_stages;
tstart = 0;
tfinal = dT;
y0 = [0; 0];
y_final = [200;0];

refine = 4;
options1 = odeset('Events',@events1,'OutputSel',1,'RelTol',tol,'AbsTol',tol/10,'Refine',refine);

v_1 = 10;
v_max = 25;
tout = tstart;
yout = y0.';
teout = [];
yeout = [];
ieout = [];

mode = 0;
mode_opt = [mode];
for i = 1:N_stages
    u = u_opt(i);
    % Solve until the first terminal event.
    t = tstart;
    tout_i = [];
    yout_i = [];

    while abs(t(end)-tfinal)>tol
        if ~isequal(tstart,tfinal)
            if mode == 0
                [t,y,te,ye,ie] = ode45(@(t,y) [y(2); u ],[tstart tfinal],y0,options1);
                teout = [teout; te];          % Events at tstart are never reported.
                yeout = [yeout; ye];
            else
                if hybrid_dynamics
                  [t,y,te,ye,ie] = ode45(@(t,y) [y(2); 3*u ],[tstart tfinal],y0,options1);
                else
                   [t,y,te,ye,ie] = ode45(@(t,y) [y(2); u ],[tstart tfinal],y0,options1);
                end
                teout = [teout; te];          % Events at tstart are never reported.
                yeout = [yeout; ye];
            end
            tout_i = [tout_i ;t];
            yout_i = [yout_i ;y];
        end

        if isequal(ie,1)
            if abs(te-tstart)<tol
%                 tstart = tfinal;
                tstart = te+2*tol;
                y0 = y(end,:);
            else
                mode = 1-mode;
                tstart = te;
                y0 = y(end,:);
            end
            ie = [];
            ye = nan;
        end
    end

    mode_opt  = [mode_opt;mode];
    tout = [tout; tout_i];
    yout = [yout; yout_i];
    ieout = [ieout; ie];

    % Set the new initial conditions
    y0 = y(end,:);
    tstart = t(end);
    tfinal = tstart + dT;
end
if 1
figure
plot(tout,yout(:,2),'LineWidth',1.5)
hold on
yline(v_max,'r--','LineWidth',1.5)
yline(v_1,'k--','LineWidth',1.5)
yline(0,'k-','LineWidth',1.5)
xlabel('$t$','Interpreter','latex')
ylabel('$v(t)$','Interpreter','latex')
ylim([min(yout(:,2))-1 max(yout(:,2))+1])
xlim([0 T_opt])
%     grid on
try
    plot(teout,yeout(:,2),'rx')
catch
end
xlabel('time');
ylabel('velocity');
for ii = 1:N_stages
    xline(dT*ii,'k:');
end
end
% mode_opt
error = norm(y_final-yout(end,1:2)');
fprintf('Terminal error is : %2.2e \n',error);
% --------------------------------------------------------------------------

    function [value,isterminal,direction] = events1(t,y)
        % Locate the time when height passes through zero in a decreasing direction
        % and stop integration.
        value = y(2)-10;     % detect height = 0
        isterminal = 1;   % stop the integration
        direction = 0;   % postive direction
    end
% --------------------------------------------------------------------------

end


