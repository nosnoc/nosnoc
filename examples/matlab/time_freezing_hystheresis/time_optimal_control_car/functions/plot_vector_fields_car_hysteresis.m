%% ilustration of the hystheresis time freezing



scale_factor = 3;
x1 = v1/scale_factor;
x2 = v2/scale_factor;
u = sign(u_max)*0.1;

x_prolonge = 2;
%
w_min = a1;
w_max = a2;



y_max = 2.4;
y_max = max(x1,x2)+4;
dx = 0.25; % resolution for vector field


my_color1 = [0 0 0]; % black
my_color2 = [0.8500 0.3250 0.0980]; % red
my_color3 = [0 0.4470 0.7410]; 	 % blue

my_color4 = [0.9290 0.6940 0.1250]; % orange
my_color5 = [0.4660 0.6740 0.1880]; % green
my_color6 = [0.6350 0.0780 0.1840]; % blood red

my_color6 = my_color3;

boundary_linewidth = 1;
charahteristic_linewidth = 2.0;
charahteristic_marker_size = 4;
plot_all_details = 1;
plot_quvier = 1;

figure

% plot caracteristic
% tt1 = x1:0.1:x2+x_prolonge ;
% plot(tt1,tt1*0+w_max,color=my_color6 ,linewidth=charahteristic_linewidth);
% hold on
% tt1 = x1-x_prolonge:0.1:x2;
% plot(tt1,tt1*0+w_min,color=my_color6 ,linewidth=charahteristic_linewidth);
% plot(x2,w_min,'o','Color',my_color6,'MarkerEdgeColor',my_color6,'MarkerFaceColor',my_color6,MarkerSize=charahteristic_marker_size);
% plot(x1,w_max,'o','Color',my_color6,'MarkerEdgeColor',my_color6,'MarkerFaceColor',my_color6,MarkerSize=charahteristic_marker_size);
% tt1 = w_min:0.1:w_max;
% plot(tt1*0+x1,tt1,color=my_color6 ,linewidth=charahteristic_linewidth,LineStyle=":");
% plot(tt1*0+x2,tt1,color=my_color6 ,linewidth=charahteristic_linewidth,LineStyle=":");

plot(x2_opt/scale_factor,x4_opt,LineWidth= 3)
hold on

    % plot region boundaries
    tt1 = -y_max:0.1:y_max;
    % horizontal boundaries
    plot(tt1,tt1*0+w_min,color=my_color1,linewidth=boundary_linewidth ,LineStyle="-");
    plot(tt1,tt1*0+w_max,color=my_color1,linewidth=boundary_linewidth ,LineStyle="-");

    % vecrtical bounds
    plot(tt1*0+x1,tt1,color=my_color1,linewidth=boundary_linewidth ,LineStyle="--");
    plot(tt1*0+x2,tt1,color=my_color1,linewidth=boundary_linewidth ,LineStyle="--");


    % diagonal
    tt1 = -y_max:0.1:y_max;
    k_lin = (w_min-w_max)/(x2-x1);
    a_lin = w_min-k_lin*x2;
    plot(tt1,k_lin*tt1+a_lin,color=my_color1,linewidth=boundary_linewidth ,LineStyle="-");


    transparency_factor = 0.3;
    % lower left - lower right - upper right - upper left
    % Region 1: DAE forming for w = 0;
    x = [-y_max y_max x2 -y_max];
    y = [k_lin*y_max+a_lin k_lin*y_max+a_lin w_min w_min];
    patch(x,y,my_color5,'FaceAlpha',transparency_factor)

    % Region 2 and Region 3 ( push down dynamics) towards w = 0
    x = [-y_max x2 x1 -y_max];
    y = [w_min w_min w_max w_max];
    patch(x,y,my_color2,'FaceAlpha',transparency_factor)
    x = [-y_max x1 -y_max];
    y = [w_max w_max k_lin*(-y_max)+a_lin];
    patch(x,y,my_color2,'FaceAlpha',transparency_factor-0.1)

    % Region 4

    x = [x1 y_max y_max -y_max];
    y = [w_max w_max k_lin*(-y_max)+a_lin k_lin*(-y_max)+a_lin];
    patch(x,y,my_color5,'FaceAlpha',transparency_factor)

    % Region 5
    x = [x2 y_max y_max x1];
    y = [w_min w_min w_max w_max];
    patch(x,y,my_color4,'FaceAlpha',transparency_factor)

    %     Region 6
    x = [x2 y_max y_max ];
    y = [w_min w_min k_lin*(y_max)+a_lin];
    patch(x,y,my_color4,'FaceAlpha',transparency_factor-0.1)


% u = -0.1;
a_push = 0.05;
    % plot vector fields
    if plot_quvier
        % region r1
        
        [X,Y] = meshgrid(-y_max:dx:y_max,-y_max:dx:y_max);
        sigma = 1e-5;
        M1 = 0.5*(1+tanh((Y-k_lin*X-a_lin)./sigma));
        M2 = 0.5.*(1+tanh((Y-w_min)./sigma));
        M3 = 0.5.*(1+tanh((Y-w_max)./sigma));
        indicator_R1 = (1-M1).*(1-M2);
        U = indicator_R1.*(u);

        if model.linear_auxiliary_dynamics
         V = -a_push*indicator_R1.*(X-x2);
        else
            V = a_push*indicator_R1;
        end
        quiver(X,Y,U,V,0,Color=my_color5);

        % regions r2 and r3
        indicator_R23 = (1-M1).*(M2);
        U = indicator_R23.*0;
        if model.linear_auxiliary_dynamics
            V = a_push*indicator_R23.*(X-x2);
        else
            V = -a_push*indicator_R23;
        end
        quiver(X,Y,U,V,0,Color=my_color2)

        % regions r4
        indicator_R4 = (M1).*(M3);
        U = indicator_R4.*(2*u);
        if model.linear_auxiliary_dynamics
            V = -a_push*indicator_R4.*(X-x1);
        else
            V = -a_push*indicator_R4;
        end
        
        quiver(X,Y,U,V,0,Color=my_color5)

        % regions r56
        indicator_R56 = (M1).*(1-M3);
        U = indicator_R56.*0;
        if model.linear_auxiliary_dynamics
            V = a_push*indicator_R56.*(X-x1);
        else
            V = a_push*indicator_R56;
        end
        quiver(X,Y,U,V,0,Color=my_color4)
    end
%         text(-1,-0.25,'$R_1$',Interpreter='latex',FontSize=15)
%         text(-2,0.64,'$R_2$',Interpreter='latex',FontSize=15)
%         text(-1.5,1.8,'$R_3$',Interpreter='latex',FontSize=15)
%         text(0.0,1.8,'$R_4$',Interpreter='latex',FontSize=15)
%         text(1.5,1,'$R_5$',Interpreter='latex',FontSize=15)
%         text(1.5,-0.25,'$R_6$',Interpreter='latex',FontSize=15)


xlim([-y_max y_max])
ylim([w_min-0.5 w_max+0.5])


xlabel('$x$','Interpreter','latex')
ylabel('$w$','Interpreter','latex')
% xticklabels([])
% yticklabels([])

%% vector fields in x
function y = f(x)
%     y = 0.14*cos(x);
y = u;
end

function y = g1(x,y)
% for w = 0
%    y = 0.15*cos(x);
% y = -0.01*x+0.03+(1-y)*0.2;
y = 1*0;
end
function y = g2(x,y)
% for w = 1
y = 2*0;
end



