%
%    This file is part of NOS-NOC.
%
%    NOS-NOC -- A software for NOnSmooth Numerical Optimal Control.
%    Copyright (C) 2022 Armin Nurkanovic, Moritz Diehl (ALU Freiburg).
%
%    NOS-NOC is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 3 of the License, or (at your option) any later version.
%
%    NOS-NOC is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with NOS-NOC; if not, write to the Free Software Foundation,
%    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%
close all
%% read data
unfold_struct(results,'base');
unfold_struct(stats,'base');
unfold_struct(model,'base');

tgrid = cumsum([0; h_vec(1:1:end)]);
x1_opt = x_res(1,:);
x2_opt = x_res(2,:);
x3_opt = x_res(3,:);
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
plot(x3_opt,x2_opt,'linewidth',1.5)
ylabel('$w(t)$','Interpreter','latex')
grid on
xlim([0 max(x3_opt)]);
subplot(326)
plot(x3_opt,x3_opt,'linewidth',1.5)
xlim([0 max(x3_opt)]);
hold on
ylim([0 max(x3_opt)])
ylabel('$t(\tau)$ ','Interpreter','latex')
xlabel('$t$ [phyisical time]','Interpreter','latex')
grid on
