 % plot scenario 
x_target = 3;
xmin = -2;
xmax = x_target+1;
temp = xmin:1:xmax;
n_holes = 3;
xc_vec = [0.5 1.5 2.5];
zc_vec = [0 0 0];
height_vec = [0.1 0.1 0.1];
width_vec = [0.5 0.5 0.5];
a_vec = width_vec/2;
b_vec = height_vec;
 area(temp ,0*temp-1,'FaceColor','k','FaceAlpha',0.2,'EdgeColor',[1 1 1]);
 hold on
for kk= 1:n_holes
                x_nodes = [xc_vec(kk)-a_vec(kk)  xc_vec(kk)+a_vec(kk) xc_vec(kk)+a_vec(kk) xc_vec(kk)-a_vec(kk)];
                y_nodes = [-5 -5 1e-2 1e-2];
                patch(x_nodes,y_nodes,'white','EdgeColor','white');
end
xline(0,'r--')
xline(x_target,'r--')
axis equal
xlim([-0.5 xmax-0.5]);
ylim([-0.5 1]);