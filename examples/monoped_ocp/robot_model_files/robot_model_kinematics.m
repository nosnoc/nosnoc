import casadi.*
 % positions
% Hip
p_hipCOM = [q(1);q(2)];
% CoM
p_ThighCOM = p_hipCOM + sThigh*[-sin(q(3));-cos(q(3))]; 
p_ShankCOM = p_hipCOM + lThigh* [-sin(q(3));-cos(q(3))]+ sShank*[-sin(q(3)+q(4));-cos(q(3)+q(4))];
% Interesting Points
% p_head = p_hipCOM + lHead*[-sin(q(3));-cos(q(3))];
base_angle = pi/4*0;
p_head = p_hipCOM + lHead*[sin(base_angle);cos(base_angle)];
p_knee = p_hipCOM + lThigh* [-sin(q(3));-cos(q(3))]; 
% p_knee = p_hipCOM + lThigh* [-sin(-q(3));-cos(-q(3))]; 
p_foot = p_hipCOM + lThigh* [-sin(q(3));-cos(q(3))]+ lShank*[-sin(q(3)+q(4));-cos(q(3)+q(4))]; 

% velocity
v_hipCOM = p_hipCOM.jacobian(q)*v;
v_ThighCOM = jacobian(p_ThighCOM,q)*v;
v_ShankCOM = jacobian(p_ShankCOM,q)*v;
v_foot = jacobian(p_foot,q)*v;
v_knee = jacobian(p_knee,q)*v;
% casadi functins
p_head_fun = Function('p_head_fun',{x},{p_head});
p_foot_fun = Function('p_foot_fun',{x},{p_foot});
p_knee_fun = Function('p_knee_fun',{x},{p_knee});
v_foot_fun = Function('v_foot_fun',{x},{v_foot});
v_knee_fun = Function('v_knee_fun',{x},{v_knee});
v_hipCOM_fun  = Function('v_hipCOM_fun',{x},{v_hipCOM});


% intertia matrix
try
    M_fun = Function.load('M_fun');
    h_forces_fun= Function.load('h_forces_fun');
catch
    M_fun = Function.load('robot_model_files/M_fun');
    h_forces_fun= Function.load('robot_model_files/h_forces_fun');
end
M = M_fun([x]);
h_forces = h_forces_fun([x]);