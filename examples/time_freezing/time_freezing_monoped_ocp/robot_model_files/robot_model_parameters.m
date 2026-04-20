%
mHip = 3.975; % mass of hip
mThigh = 1.782; %
mShank = 0.548;
%lengts
lBH = 0.043;  % distance between base and hip joint
lThigh = 0.2;
lShank = 0.2;
lHead = 0.05;
% center of masses distances
sBM = 0.02;  % distance between base and CoG of main body
sThigh = 0.016; % distance between hip joint and CoG of thigh *can be negative if z is positive?
sShank = 0.1;  % distance between knee joint and CoG of shank
rf = 0.028; % radius of foot
IyThugh = 0.001; % kgm2 inertia of thigh w.r.t. CoG about z-axis
IyShank = 0.0032;% inertia of shank w.r.t. CoG about z-axis
g = 9.81;