clear all 
clc
import casadi.*
%% Integrator settings
omega = 2*pi;
A1 = [1 omega;...
    -omega 1];
A2 = [1 -omega;...
    omega 1];

%%
h =  0.004;

integrator = 'idas';
tol = 1e-9;

% integrator = 'collocation';
% opts = struct('tf', h, 'reltol',tol*1e-1,'abstol',tol);

collocation_scheme = 'radau';
interpolation_order = 1;
number_of_finite_elements = 1;

% interpolation_order = 2;
% number_of_finite_elements = 1;

%%  load model (variable delcaration, parameters ect.)
x0 = [exp(-1);0];

x1 = MX.sym('x1');
x2 = MX.sym('x2');
x = [x1;x2];
n_x = 2;
f_x = A2*x;

psi = x1^2+x2^2-1;
sigma = 1e-6;
alpha = (1+tanh(psi/sigma))/2;
f_x = alpha*A1*x+(1-alpha)*A2*x;

%% create casadi integrators

% integrator struct
dae = struct('x', x, 'ode', f_x);
T = 2; 
N = T/h;

opts = struct('tf', h, 'reltol',tol*1e-1,'abstol',tol);

if isequal(integrator,'collocation')
%     collocation_scheme = 'radau';
%     interpolation_order = 1;
%     number_of_finite_elements = 20;
    
%     interpolation_order = 2;
%     number_of_finite_elements = 2;
    
    opts = struct('tf', h,'collocation_scheme',collocation_scheme,'interpolation_order',interpolation_order,...
           'number_of_finite_elements',number_of_finite_elements);
end




%% simulation loop
% % Integrator
I = casadi.integrator('I', integrator, dae, opts);

% %%
x_res = [x0];
tic
for ii = 1:N
    try
        res = I('x0', x0);       
    catch
        fprintf('Failed at ii = %d, t  = %5.2f \n',ii,ii*h);
        error('divergence');
    end
    xf = full(res.xf);
    %     
    x_res = [x_res,xf];
    x0 = xf;
end
t_cpu_sim = toc;
fprintf('total simulation time %3.2f s \n',t_cpu_sim);
time = 0:h:T;
%% read 
x1_opt = x_res(1,:);
x2_opt = x_res(2,:);

figure
subplot(121)
% fimplicit(@(x,y) (x-0).^2+(y-0).^2-1^2, [-3 3],'linewidth',2.0,'color',0.55*ones(3,1))
theta = 0:0.1:2*pi;
x = sin(theta);
y = cos(theta);
plot(x1_opt,x2_opt,'k','linewidt',1.0);
hold on
plot(x,y,'linewidth',1.0,'color',0.4*ones(3,1))
hold on
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
axis equal
grid on
% [X,Y] = meshgrid(-3:0.35:3,-3:0.35:3);
[X,Y] = meshgrid(-3:0.5:3,-3:0.5:3);
[U,V] = oscilator_eval(X,Y);
quiver(X,Y,U,V,'Color',0.75*ones(3,1),'LineWidth',0.8,'AutoScaleFactor',3);
% [X,Y] = meshgrid(-1:0.1:1,-1:0.1:1);
% [U,V] = oscilator_eval(X,Y);
% quiver(X,Y,U,V,'Color',0.6*ones(3,1),'LineWidth',1.2,'AutoScaleFactor',2);
xlim([-exp(1) exp(1)]*1.01)
ylim([-exp(1) exp(1)]*0.9)
%
subplot(122)
% figure
plot(time,x1_opt,'color',0.6*ones(3,1),'linewidth',1.0);
hold on
plot(time,x2_opt,'color',0*ones(3,1),'linewidth',1.0);
xlabel('$t$','interpreter','latex');
ylabel('$x(t)$','interpreter','latex');
legend({'$x_1(t)$','$x_2(t)$'},'interpreter','latex');
grid on
% axis equal
%% Plots 
norm([x1_opt(end);x2_opt(end)]-[1;0])

