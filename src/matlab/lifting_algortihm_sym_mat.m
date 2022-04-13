clc
clear all
n_f = 6;
n_c = 5;

% it works
S = [1 0 0 0 0;...
    -1 1 0 0 0;...
    -1 -1  1 0 0;...
    -1 -1 -1 1 0;...
    -1 -1 -1 -1 1;...
    -1 -1 -1 -1 -1;...
    ];


% S = [1 0;...
%     -1 1;...
%     -1 -1];
%
%
% S = [1 0 0;...
%     -1 1 0;...
%     -1 -1 -1];
%
% %it fails
% % %
% S = [1 1 1;...
%      1 -1 1;...
%      1 -1 -1;...
%      -1 1 1;...
%      -1 -1 1;...
%      -1 -1 -1;...
%      ];
% % % %
% S = [1 1;...
%      1 -1;...
%      -1 1;...
%      -1 -1];
%


[n_f,n_c] = size(S);
n_R = sum(abs(S),2);

alpha = sym('alpha', [n_c 1]);
theta = sym('theta', [n_f 1]);

% alpha = MX.sym('alpha',n_c);
% theta = MX.sym('theta',n_f);
%%
s = 0;
ii = 1;
chi = (1-s)/2+s*alpha(ii);

%%
% pseudo code
% are beta needed?
% yes, produce the beta
% no ,
% produce temp
% add temp and remove what is too much
g_lift_beta = [];
g_lift_theta = [];
beta = [];

temp = ones(n_f,1);
n_depth = 2;
temp_S  = S;
ind_progress = 0;
for ii = 1:n_c
    n_R = sum(abs(temp_S),2);
    ind_done = find(n_R >= ii);
    temp = (temp).*((1-temp_S(ind_done,ii))/2+temp_S(ind_done,ii).*alpha(ii));
    ind_done= find(n_R == ii);
    if ~isempty(ind_done)
        g_lift_theta = [g_lift_theta ;theta(ind_progress+ind_done)-(temp(ind_done))];
        temp(ind_done) = [];
        temp_S(ind_done,:) = [];
        ind_progress = ind_progress+ind_done(end);
    end
    % start defining betas;
    if ii >= n_depth && ii< n_c
        [temp_S_red,temp_S_IA,temp_S_IC] = unique(temp_S(:,ii) ,'rows');
        n_beta_ii = size(temp_S_red,1);
        beta_temp = sym(['beta_' num2str(ii+1-n_depth)], [n_beta_ii 1]);
        beta = [beta;beta_temp]
        g_lift_beta = [g_lift_beta;beta_temp - temp(temp_S_IA)];
        temp = beta_temp(temp_S_IC);
    end
end
beta
g_lift_beta
g_lift_theta
