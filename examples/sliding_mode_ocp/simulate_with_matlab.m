% model equations
u_opt = [0.0287214975870860	-3.25252098226141e-07	-9.91402433240917e-08	1.86445793774576e-05	-0.226932588887457	-2.17832668754677
    0.00190634490399469	0.00704916496394133	-0.0453974298042458	0.169645233740009	-1.80539211250271	0.0654278258088038];
% Switching Functions
p = 2; a = -0.1; a1 = 0;
b = -0.05; q = 3;
%
tspan = [0 4/6];
y0 = [2*pi/3;pi/3;0;0];
sigma_int = 1e-12;
tol = sigma_int/2;
options = odeset('RelTol', tol, 'AbsTol', tol/10);
y_integrator = [];
t_integrator = [];
for ii = 1:6

        [t,y_res] = ode15s(@(t,y) ...
            ((0.5*(1+tanh((y(1)+a*(y(2)-a1)^p)/sigma_int))))*([-1+y(3);0;u_opt(1,ii);u_opt(2,ii)])+(1-(0.5*(1+tanh((y(1)+a*(y(2)-a1)^p)/sigma_int))))*([1+y(3);0;u_opt(1,ii);u_opt(2,ii)])+...
            ((0.5*(1+tanh((y(2)+b*y(1)^q)/sigma_int))))*([0;-1+y(4);u_opt(1,ii);u_opt(2,ii)])+(1-(0.5*(1+tanh((y(2)+b*y(1)^q)/sigma_int))))*([0;1+y(4);u_opt(1,ii);u_opt(2,ii)])...
            ,tspan, y0,options);
    
%         [t,y_res] = ode15s(@(t,y) ...
%             (1-(0.5*(1+tanh((y(1)+a*(y(2)-a1)^p)/sigma))))*([-1+y(3);0;u_opt(1,ii);u_opt(2,ii)])+((0.5*(1+tanh((x1+a*(x2-a1)^p)/sigma))))*([1+y(3);0;u_opt(1,ii);u_opt(2,ii)])+...
%             (1-(0.5*(1+tanh((y(2)+b*y(1)^q)/sigma))))*([0;-1+y(4);u_opt(1,ii);u_opt(2,ii)])+((0.5*(1+tanh((y(2)+b*y(1)^q)/sigma))))*([0;1+y(4);u_opt(1,ii);u_opt(2,ii)])...
%             ,tspan, y0);
    
    y0 = y_res(end,:)';
    y_integrator = [y_integrator,y_res'];
    t_integrator = [t_integrator,t'+(ii-1)*4/6];
end
%
figure
plot(y_integrator(1,:),y_integrator(2,:),'LineWidth',2)
if illustrate_regions
    hold on
    p = 2; a = -0.1; a1 = 0; 
    b = -0.05; q = 3;
    t2 = -5:0.01:5;
    plot(-a*(t2-a1).^p,t2,'k')
    hold on
    t1 = t2;
    plot(t1,-b*t1.^q,'k')
    grid on
    axis equal
    xlim([-1.5 1.5])
    ylim([-1.5 1.5])
end
norm(x_target-y_integrator(1:2,end))
norm(x_res_optimizer(1:2,end)-y_integrator(1:2,end))
norm(x_res_integrator(1:2,end)-y_integrator(1:2,end))
norm(x_res_integrator(1:2,end)-x_target)
