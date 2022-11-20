%% sliding mode ocp  - plots for paper
unfold_struct(results,'base');
% unfold_struct(model,'base');
unfold_struct(settings,'base');
%% geometric trajectory and time trajectorz
figure
subplot(121)
plot(x_res_integrator(1,:),x_res_integrator(2,:),'LineWidth',2)
hold on
grid on
if illustrate_regions
    hold on
    p = 2; a = 0.15; a1 = 0;
    b = -0.05; q = 3;
    t2 = -5:0.01:5;
    plot(-a*(t2-a1).^p,t2,'k')
    hold on
    t1 = t2;
    plot(t1,-b*t1.^q,'k')
    grid on
%     axis equal
    xlim([-2.2 2.2])
    ylim([-2.2 2.2])
end
plot(x_target(1),x_target(2),'rx')
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
legend({'$x(t)$','$\varphi_1(x)=0$','$\varphi_2(x)=0$'},'Interpreter','latex','Location','best');
subplot(122)
plot(t_grid_integrator,x_res_integrator(1:2,:))
grid on
xlabel('$t$','interpreter','latex');
ylabel('$x(t)$','interpreter','latex');
legend({'$x_1(t)$','$x_2(t)$'},'Interpreter','latex');
%% control functions
u1_opt = u_opt(1,:);
u2_opt = u_opt(2,:);
u_bar = max(abs(u_opt(:)));
v_bar = max(max(abs(x_res_optimizer(3:4,:))));
t_grid_control = t_grid_optimizer(1:N_finite_elements:end);
figure
subplot(121)
plot(t_grid_optimizer,x_res_optimizer(3,:),'LineWidth',1.5)
hold on
plot(t_grid_optimizer,x_res_optimizer(4,:),'LineWidth',1.5)
xlabel('$t$','interpreter','latex');
ylabel('$v(t)$','interpreter','latex');
ylim([-1.1*v_bar 1.1*v_bar])
for ii = 1:6
    xx = -10:2:10;
    plot(xx*0+t_grid_control(ii+1),xx,'k--')
end
grid on
legend({'$v_1(t)$','$v_2(t)$'},'Interpreter','latex','Location','best');
subplot(122)
stairs(t_grid_control,[u1_opt,nan],'LineWidth',1.5)
hold on
stairs(t_grid_control,[u2_opt,nan],'LineWidth',1.5)
for ii = 1:6
    xx = -10:2:10;
    plot(xx*0+t_grid_control(ii+1),xx,'k--')
end
ylim([-1.1*u_bar 1.1*u_bar])
grid on
xlabel('$t$','interpreter','latex');
ylabel('${u}(t)$','interpreter','latex');
legend({'${u}_1(t)$','${u}_2(t)$'},'Interpreter','latex','Location','best');

%% indicator functio
if 0
if step_equilibration
    nu_opt = full(model.nu_fun(results.w_opt));
else
    nu_opt = zeros(N_finite_elements*6-1,1);
end
nu_opt_scaled = tanh(nu_opt/0.001);
figure
subplot(121)
semilogy(nu_opt)
grid on
subplot(122)
% figure
% stairs(nu_opt_scaled)
plot(nu_opt_scaled)
ylim([-0.05 1.05])
grid on
end
%%
if 0
figure
subplot(131)
try

semilogy(1:N_finite_elements*6,[nan;nu_opt],'LineWidth',1.5);
hold on
% for ii = 1:N_finite_elements:N_finite_elements*6
%     xx = logspace(1e-16,1e1,10);
%     hold on
%     semilogy(xx*0+ii,xx,'k:')
% end
catch
    semilogy(nu_opt)
end
xlabel('Grid point','interpreter','latex');
ylabel('$\eta_{n,k}(\cdot)$','interpreter','latex');
xlim([1 N_finite_elements*6])
grid on

%
h_opt =results.h_opt;
h_bar = model.h*2;
h_bar = max(results.h_opt);
h_nominal = model.h/N_finite_elements;
% figure
subplot(132)
stairs(h_opt,'LineWidth',1.5)
% grid on
hold on
for ii = 1:N_finite_elements:N_finite_elements*6
    xx = 0:1:2;
    plot(xx*0+ii,xx,'k:')
end
xx = 0:2:N_finite_elements*6;
plot(xx,xx*0,'k--')
plot(xx,xx*0+h_nominal *2,'k--')
plot(xx,xx*0+h_nominal ,'k--')
xlim([1 N_finite_elements*6])
ylim([-0.01 h_nominal*2.1])
xlabel('Finite element','interpreter','latex');
ylabel('$h_{k,n}$','interpreter','latex');
subplot(133)
stairs(h_opt,'LineWidth',1.5)
% grid on
hold on
for ii = 1:N_finite_elements:N_finite_elements*6
    xx = 0:1:2;
    plot(xx*0+ii,xx,'k:')
end
xx = 0:2:N_finite_elements*6;
plot(xx,xx*0,'k--')
plot(xx,xx*0+h_nominal *2,'k--')
plot(xx,xx*0+h_nominal ,'k--')
xlim([1 N_finite_elements*6])
ylim([-0.01 h_nominal*2.1])
xlabel('Finite element','interpreter','latex');
% ylabel('$h_{k,n}$','interpreter','latex');

end