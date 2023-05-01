function [Quadruped_X, Quadruped_Y] = getQuadrupedConfiguration(model, q)
x = q(1);
z = q(2);
xita_torso = q(3);
xita_thigh_1 = q(4);
xita_calf_1 = q(5);
xita_thigh_3 = q(6);
xita_calf_3 = q(7);


l_torso = model.linkLength(1);
l_thigh = model.linkLength(2);
l_calf = model.linkLength(3);

% torso
torso_X = [x;...
    x + l_torso * sin(xita_torso)];
torso_Y = [z;...
    z - l_torso * cos(xita_torso)];

% leg 1
leg_1_X = [x; ...
    x + l_thigh * sin(xita_thigh_1);...
    x + l_thigh * sin(xita_thigh_1) + l_calf * sin(xita_calf_1)];
leg_1_Y = [z;...
    z - l_thigh * cos(xita_thigh_1);...
    z - l_thigh * cos(xita_thigh_1) - l_calf * cos(xita_calf_1)];


% leg 3
leg_3_X = [x + l_torso * sin(xita_torso);...
    x + l_torso * sin(xita_torso) + l_thigh * sin(xita_thigh_3);...
    x + l_torso * sin(xita_torso) + l_thigh * sin(xita_thigh_3) + l_calf * sin(xita_calf_3)];
leg_3_Y = [z - l_torso * cos(xita_torso);...
    z - l_torso * cos(xita_torso) - l_thigh * cos(xita_thigh_3);...
    z - l_torso * cos(xita_torso) - l_thigh * cos(xita_thigh_3) - l_calf * cos(xita_calf_3)];

Quadruped_X = [torso_X; leg_1_X; leg_3_X];
Quadruped_Y = [torso_Y; leg_1_Y; leg_3_Y];
end
