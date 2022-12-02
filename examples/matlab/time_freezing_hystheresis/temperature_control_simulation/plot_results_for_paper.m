% Copyright 2022 Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% This file is part of NOSNOC.

% The 2-Clause BSD License

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
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

saveas(gcf,'numerical_vs_physical_time')

%% complementarity stats
figure
semilogy(complementarity_stats+1e-16,'k',LineWidth=1.5)
grid on
xlabel('$\tau$ [numerical time]','Interpreter','latex')
ylabel('Complementarity residual','Interpreter','latex')