%%
rho = -2:0.001:10;
y1 = -exp(-rho);
y2 = exp(-rho);
figure
plot(y1,rho);
xline(0)
hold on
grid on
plot(y2,rho);
xlabel('\sigma')
ylabel('\rho')


%%
figure
plot(rho,y1);
xline(0)
hold on
grid on
plot(rho,y2);
ylabel('\sigma')
xlabel('\rho')
%%
figure
plot(-log(exp(-rho)))


%%
rho = -10:0.001:10;
y1 = -exp(rho);
y2 = exp(rho);
figure
plot(y1,rho);
xline(0)
hold on
grid on
plot(y2,rho);
xlabel('\sigma')
ylabel('\rho')

%%
mu = 1e-16
rho = 0.01:0.01:10;
plot(rho,-(1/mu)*log(rho))

mu = 1e-2;
rho = 0.01:0.01:10;
plot(rho,-(1/mu)*log(rho))

%%
lambda = 1;
rho = linspace(-1,10,200);
sigma = linspace(-5,5,200);
plot(exp(-lambda*rho),rho,'k')
xline(0,'k')
hold on
grid on
[RHO,SIGMA] = meshgrid(rho,sigma);
Z = -RHO.^2;
% contour(SIGMA,RHO,Z,30);
figure
t = 100;
for i = 1:20
    t = t/2;
    clf
    plot(exp(-lambda*rho),rho,'k')
    hold on
    xline(0,'k')
    grid on
    Z = -RHO.^2-(1/t)*(log(max(1e-4,SIGMA))+log(max(1e-4,exp(-RHO)-SIGMA)));
    contour(SIGMA,RHO,Z)
    colorbar
    caxis([-1e2 1e3])

%     set(gca,'ColorScale','log')
    xlim([-1 5])
    ylim([-1 12])
    ylabel('\rho')
    xlabel('\sigma')
    pause (1)

%     sigma = t;
%     rho = 10;
%     obj = -rho.^2-(1/t)*(log(max(1e-4,sigma))+log(max(1e-4,exp(-rho)-sigma)))
end




