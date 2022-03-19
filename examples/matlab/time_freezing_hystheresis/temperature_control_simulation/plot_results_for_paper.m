close all
unfold_struct(results,'base');
tgrid = cumsum([0; h_vec(1:1:end)]);
% tgrid = t;
x1_opt = x_res(1,:);
x2_opt = x_res(2,:);
x3_opt = x_res(3,:);
% tgrid = tgrid();
y1 = 18;
y2 = 20;
theta1 = 18;
theta0 = 15;
theta2 = 20;

%% in numerical time
t_phy = [x3_opt,nan];
dh = diff([tgrid;nan]');
dt_phy = diff(t_phy)./dh;
fixer = nan;
dt_phy(dt_phy==-inf) = fixer ;
dt_phy(dt_phy==inf) = fixer ;
ddt_phy = diff([dt_phy,nan]);

t_criterion = max(abs(dt_phy))/2;
dt_criterion = max(abs(ddt_phy))/2;
ind_t = find(abs(dt_phy)>t_criterion);
ind_t_complement = find(abs(dt_phy)<=t_criterion);
ind_dt = find(abs(ddt_phy)>dt_criterion);

time_frozen = tgrid*0;
time_frozen(ind_t_complement)=1;

break_points = tgrid(ind_dt);

%% Split in two plots
figure
subplot(311)
plot(tgrid,x1_opt,'linewidth',1.5)
hold on
for ii = 1:length(break_points)/2
    xx = [break_points(2*ii-1) break_points(2*ii) break_points(2*ii) break_points(2*ii-1)];
    yy = [0 0 max(x1_opt)+3 max(x1_opt)+3];
    patch(xx,yy,'red','facealpha',0.08)
end
hold on
plot(tgrid,tgrid*0+theta1,'k--')
plot(tgrid,tgrid*0+theta2,'k--')
ylim([theta0-1,theta2+1])
xlim([0 T_sim])
% xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('$x(\tau)$','Interpreter','latex')
grid on 
subplot(312)
plot(tgrid,x2_opt,'linewidth',1.5)
hold on
% area(tgrid,time_frozen*(max(x2_opt)+3),FaceAlpha=0.1)
for ii = 1:length(break_points)/2
    xx = [break_points(2*ii-1) break_points(2*ii) break_points(2*ii) break_points(2*ii-1)];
    yy = [0 0 max(x2_opt)+3 max(x2_opt)+3];
    patch(xx,yy,'red','facealpha',0.08)
end
ylim([0 1.2])
% xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('$w(\tau)$','Interpreter','latex')
grid on
xlim([0 T_sim])
subplot(313)
plot(tgrid,x3_opt,'linewidth',1.5)
hold on
% area(tgrid,time_frozen*(max(x3_opt)+1),FaceAlpha=0.1)
for ii = 1:length(break_points)/2
    xx = [break_points(2*ii-1) break_points(2*ii) break_points(2*ii) break_points(2*ii-1)];
    yy = [0 0 max(x3_opt)+3 max(x3_opt)+3];
    patch(xx,yy,'red','facealpha',0.08)
end
ylim([0 max(x3_opt)])
xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('$t(\tau)$ ','Interpreter','latex')
grid on
xlim([0 T])
%% Phyisical time
figure
subplot(211)
plot(x3_opt,x1_opt,'linewidth',1.5)
hold on
plot(x3_opt,x3_opt*0+theta1,'k--')
plot(x3_opt,x3_opt*0+theta2,'k--')
ylim([theta0-1,theta2+1])
ylabel('$x(t)$','Interpreter','latex')
xlabel('$t$ [phyisical time]','Interpreter','latex')
grid on
xlim([0 max(x3_opt)]);
subplot(212)
% plot(x3_opt,1-x2_opt,'linewidth',1.5)
plot(x3_opt,x2_opt,'linewidth',1.5)
% xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('$w(t)$','Interpreter','latex')
xlabel('$t$ [phyisical time]','Interpreter','latex')
grid on
xlim([0 max(x3_opt)]);

%% all toghether
figure
subplot(321)
plot(tgrid,x1_opt,'linewidth',1.5)
hold on
for ii = 1:length(break_points)/2
    xx = [break_points(2*ii-1) break_points(2*ii) break_points(2*ii) break_points(2*ii-1)];
    yy = [0 0 max(x1_opt)+3 max(x1_opt)+3];
    patch(xx,yy,'red','facealpha',0.08)
end
hold on
plot(tgrid,tgrid*0+theta1,'k--')
plot(tgrid,tgrid*0+theta2,'k--')
ylim([theta0-1,theta2+1])
xlim([0 T_sim])
% xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('$x(\tau)$','Interpreter','latex')
grid on 
subplot(323)
plot(tgrid,1-x2_opt,'linewidth',1.5)
hold on
% area(tgrid,time_frozen*(max(x2_opt)+3),FaceAlpha=0.1)
for ii = 1:length(break_points)/2
    xx = [break_points(2*ii-1) break_points(2*ii) break_points(2*ii) break_points(2*ii-1)];
    yy = [0 0 max(x2_opt)+3 max(x2_opt)+3];
    patch(xx,yy,'red','facealpha',0.08)
end
ylim([0 1.2])
% xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('$w(\tau)$','Interpreter','latex')
grid on
xlim([0 T_sim])
subplot(325)
plot(tgrid,x3_opt,'linewidth',1.5)
hold on
% area(tgrid,time_frozen*(max(x3_opt)+1),FaceAlpha=0.1)
for ii = 1:length(break_points)/2
    xx = [break_points(2*ii-1) break_points(2*ii) break_points(2*ii) break_points(2*ii-1)];
    yy = [0 0 max(x3_opt)+3 max(x3_opt)+3];
    patch(xx,yy,'red','facealpha',0.08)
end
ylim([0 max(x3_opt)])
xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('$t(\tau)$ ','Interpreter','latex')
grid on
xlim([0 T_sim])
subplot(322)
plot(x3_opt,x1_opt,'linewidth',1.5)
hold on
plot(x3_opt,x3_opt*0+theta1,'k--')
plot(x3_opt,x3_opt*0+theta2,'k--')
ylim([theta0-1,theta2+1])
ylabel('$x(t)$','Interpreter','latex')
grid on
xlim([0 max(x3_opt)]);
subplot(324)
% plot(x3_opt,1-x2_opt,'linewidth',1.5)
plot(x3_opt,x2_opt,'linewidth',1.5)
% xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('$w(t)$','Interpreter','latex')
grid on
xlim([0 max(x3_opt)]);
subplot(326)
plot(x3_opt,x3_opt,'linewidth',1.5)
xlim([0 max(x3_opt)]);
hold on
% area(tgrid,time_frozen,FaceAlpha=0.1)
ylim([0 max(x3_opt)])
ylabel('$t(\tau)$ ','Interpreter','latex')
xlabel('$t$ [phyisical time]','Interpreter','latex')
grid on
%% phase plot
figure
x1_phy = x1_opt;
x1_phy(ind_t_complement) = nan;
x2_phy = x2_opt;
x2_phy(ind_t_complement) = nan;
plot(x1_opt,x2_opt,'LineWidth',1.5)
hold on
plot(x1_phy,x2_phy,'LineWidth',3.0)
grid on
ylim([-0.1 1.1])
ylabel('$w$','Interpreter','latex')
xlabel('$x$ ','Interpreter','latex')
tt = 14:25;
plot(tt,tt*0,'k')
plot(tt,tt*0+1,'k')
xlim([theta0 theta2+2])
ylim([-0.2 1.2])

w0 = 0; w1 = 1;
k =  (w0-w1)/(theta2-theta1);
a = w1-k*theta1;
a = w0-k*theta2;
plot(tt,k*tt+a,'k')

text(theta1-1,-0.1,'$R_1$',Interpreter='latex',FontSize=15)
text(theta1-1,0.5,'$R_2$',Interpreter='latex',FontSize=15)
text(theta1-1,1.1,'$R_3$',Interpreter='latex',FontSize=15)

text(theta2+1,1.1,'$R_4$',Interpreter='latex',FontSize=15)
text(theta2+1,0.5,'$R_5$',Interpreter='latex',FontSize=15)
text(theta2+1,-0.1,'$R_6$',Interpreter='latex',FontSize=15)
%% plots in phyisicaltime


% figure
% subplot(211)
% plot(x3_opt,x1_opt)
% hold on
% xlabel('$t$ [phyisical time]','Interpreter','latex')
% ylabel('$\theta(t)$ [temperature]','Interpreter','latex')
% hold on
% grid on
% yline(y1,'k--')
% yline(y2,'k--')
% subplot(212)
% plot(x3_opt,x2_opt)
% grid on
