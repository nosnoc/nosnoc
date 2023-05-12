clear all;
close all
clc;
import casadi.*

%%
plot_results = 1;
[settings] = NosnocOptions();  
settings.irk_scheme = IRKSchemes.RADAU_IIA;
settings.n_s = 1;
settings.opts_casadi_nlp.ipopt.max_iter = 1e3;
settings.print_level = 2;
settings.N_homotopy = 8;
settings.use_fesd = 1;
settings.time_freezing = 1;
settings.stagewise_clock_constraint = 0;
settings.impose_terminal_phyisical_time = 0;
settings.pss_lift_step_functions = 0;
%% model prameters
n_dim_contact  = 2;
g = 9.81;
m = 1;
l = 0.5;
J = m*l^2/12;
%% Symbolic variables and model functions
qx = SX.sym('qx',1); 
qy = SX.sym('qy',1); 
theta = SX.sym('theta',1); 
vx = SX.sym('vx',1); 
vy = SX.sym('vy',1); 
omega = SX.sym('omega',1); 

q0 = [0;0.5;-pi/6*1];
v0 = [1;0;-4*0];

q = [qx;qy;theta];
v = [vx;vy;omega];

p_com = [qx;qy];

p_left = p_com-0.5*l*[cos(theta);sin(theta)];
p_right = p_com+0.5*l*[cos(theta);sin(theta)];

model.x = [q;v]; 
model.e = 0;
model.mu = 0.0;
model.a_n = g;
model.x0 = [q0;v0]; 


model.M = diag([m,m,J]);
model.f_v = [0;-g;0];
model.f_c = [p_left(2);p_right(2)];
model.n_dim_contact = n_dim_contact ;
%% Simulation settings
N_FE = 2;
T_sim = 0.8;
N_sim = 30;
model.T_sim = T_sim;
model.N_FE = N_FE;
model.N_sim = N_sim;

settings.use_previous_solution_as_initial_guess = 0;
%% Call nosnoc Integrator
[results,stats,model] = integrator_fesd(model,settings);
%% read and plot results
unfold_struct(results,'base');
qx = x_res(1,:);
qy = x_res(2,:);
theta = x_res(3,:);
vx = x_res(4,:); 
vy = x_res(5,:);
omega = x_res(6,:);
t_clock = x_res(7,:);


%% geometric trajetcorty
x_p = [-1 1 1 -1];
y_p = [-2 -2 0 0];
if plot_results
    figure('Renderer', 'painters', 'Position', [100 100 1300 1000])
    hold on
    for ii = 1:length(qx)
        xx = [qx(ii)-l/2*cos(theta(ii)),qx(ii)+l/2*cos(theta(ii))];
        yy = [qy(ii)-l/2*sin(theta(ii)),qy(ii)+l/2*sin(theta(ii))];
        plot(xx,yy,'k','LineWidth',2);
        axis equal
        p = patch(x_p,y_p,'k');
        p.FaceAlpha = 0.2; p.EdgeColor = 'none';
        xlim([-1 1]);
        ylim([-0.2 0.8]);
        grid on
        xlabel('$q_x$','interpreter','latex');
        ylabel('$q_y$','interpreter','latex');
        pause(model.h*2)
        if ii~=length(qx)
            clf
        end
    end
end

%%
figure
subplot(311)
plot(t_clock,vx,'k');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v_x$','interpreter','latex');
subplot(312)
plot(t_clock,vy,'k');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$v_y$','interpreter','latex');

subplot(313)
plot(t_clock,omega,'k');
grid on
xlabel('$t$','interpreter','latex');
ylabel('$\omega$','interpreter','latex');

