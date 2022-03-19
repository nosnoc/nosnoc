%% ilustration of the hystheresis time freezing


% dot x = f(x) + w(g_1(x)) + (1-w)g_2(x);
% dot x = f(x) + w(g_1(x)) + (1-w)g_2(x);


% region i where w is 0 or 1,
% region ii , left upper,
% region iii , left lower,
% region iv , right upper,
% region v % right lower

% boundaries
x1 = -0.75;
x2 = 0.75;
x_prolonge = 2;
%
w_min = 0;
w_max = 1.5;

y_max = 2.4;
dx = 0.3; % resolution for vector field


my_color1 = [0 0 0]; % black
my_color2 = [0.8500 0.3250 0.0980]; % red
my_color3 = [0 0.4470 0.7410]; 	 % blue

my_color4 = [0.9290 0.6940 0.1250]; % orange
my_color5 = [0.4660 0.6740 0.1880]; % green
my_color6 = [0.6350 0.0780 0.1840]; % blood red

my_color6 = my_color3;

boundary_linewidth = 1;
charahteristic_linewidth = 1.5;
charahteristic_marker_size = 4;
plot_all_details = 1;

figure

% plot caracteristic
tt1 = x1:0.1:x2+x_prolonge ;
plot(tt1,tt1*0+w_max,color=my_color6 ,linewidth=charahteristic_linewidth);
hold on
tt1 = x1-x_prolonge:0.1:x2;
plot(tt1,tt1*0+w_min,color=my_color6 ,linewidth=charahteristic_linewidth);
plot(x2,w_min,'o','Color',my_color6,'MarkerEdgeColor',my_color6,'MarkerFaceColor',my_color6,MarkerSize=charahteristic_marker_size);
plot(x1,w_max,'o','Color',my_color6,'MarkerEdgeColor',my_color6,'MarkerFaceColor',my_color6,MarkerSize=charahteristic_marker_size);
tt1 = w_min:0.1:w_max;
plot(tt1*0+x1,tt1,color=my_color6 ,linewidth=charahteristic_linewidth,LineStyle=":");
plot(tt1*0+x2,tt1,color=my_color6 ,linewidth=charahteristic_linewidth,LineStyle=":");

if plot_all_details
    % plot region boundaries
    tt1 = -y_max:0.1:y_max;
    % horizontal boundaries
    plot(tt1,tt1*0+w_min,color=my_color1,linewidth=boundary_linewidth ,LineStyle="-");
    plot(tt1,tt1*0+w_max,color=my_color1,linewidth=boundary_linewidth ,LineStyle="-");

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



    a_push = 0.15;

    if 0
        % plot vector fields
        % region upper left (ii)

        [X,Y] = meshgrid(-y_max:dx:x1,dx+w_min:dx:y_max);
        U = 0*X;
        V = 0*Y-a_push;
        quiver(X,Y,U,V,0,Color=my_color2)
        % region lower left (iii)
        [X,Y] = meshgrid(-y_max:dx:x1-dx/2,-y_max:dx:w_min-dx/2);
        % U = 0.25*cos(X);
        U = 2*f(X)+(1-Y).*g1(X)+(Y).*g2(X);
        V = 0*Y+a_push;
        quiver(X,Y,U,V,0,Color=my_color4)


        % main region (i)
        [X,Y] = meshgrid(x1+dx:dx*2:x2-dx,-y_max:dx:y_max);
        U = 2*f(X)+(1-Y).*g1(X)+(Y).*g2(X);
        V = 0*Y;
        quiver(X,Y,U,V,0,Color=my_color5)


        % region upper right (v)

        [X,Y] = meshgrid(x2:dx:y_max,dx+w_max:dx:y_max);
        U = 2*f(X)+(1-Y).*g1(X)+(Y).*g2(X);
        V = 0*Y-0.15;
        quiver(X,Y,U,V,0,Color=my_color4)
        % region lower right (iv)
        [X,Y] = meshgrid(x2:dx:y_max,-y_max:dx:w_max-dx/2);
        U = 0*X;
        V = 0*Y+a_push;
        quiver(X,Y,U,V,0,Color=my_color2)

        % regions
        text(0,0.64,'$R_1$',Interpreter='latex',FontSize=15)

        text(-2,0.64,'$R_2$',Interpreter='latex',FontSize=15)
        text(-2,-0.25,'$R_3$',Interpreter='latex',FontSize=15)

        text(1.7,0.64,'$R_4$',Interpreter='latex',FontSize=15)
        text(1.7,w_max+0.25,'$R_5$',Interpreter='latex',FontSize=15)
    end
else
    %     grid on
end




xlim([-y_max y_max])
ylim([w_min-0.5 w_max+0.5])



xlabel('$x$','Interpreter','latex')
ylabel('$w$','Interpreter','latex')
xticklabels([])
yticklabels([])

%% vector fields in x
function y = f(x)
y = 0.15*cos(x);
end

function y = g1(x)
% for w = 0
y = 0.15*cos(x);
end
function y = g2(x)
% for w = 1
y = 0.2*sin(x);
end



