%% benchmark_globals

% dream values
NS_VALUES = [1, 2, 3, 4];
NSIM_VALUES = [30, 50, 60, 80, 90, 120, 150, 200, 250, 300, 350, 450, 520, 600, 666, 900, 1200];
NFE_VALUES = [2];
IRK_SCHEME = IRKSchemes.GAUSS_LEGENDRE;

% % test values
% NS_VALUES = [2];
% NSIM_VALUES = [150];
% NFE_VALUES = [2];
% IRK_SCHEME = IRKSchemes.GAUSS_LEGENDRE;

T_sim = .3;
g = 9.81;
R = 0.2;
k = 1e4;
l = 1;
m = 1;
e = 0.8;

% x0 = [1;2;0;0];

x0 = [.3; 1.3; 0; 0];
