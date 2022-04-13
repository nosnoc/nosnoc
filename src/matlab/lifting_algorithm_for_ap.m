
import casadi.*



n_sys = 6;
n_c = 5;


S = [1 0 0 0 0;...
    -1 1 0 0 0;...
    -1 -1  1 0 0;...
    -1 -1 -1 1 0;...
    -1 -1 -1 -1 1;...
    -1 -1 -1 -1 -1;...
    ];

S = [1 1;...
     1 -1;...
     -1 1;...
     -1 -1];

S = [1 0;...
     -1 1;...
     -1 -1];


[n_sys,n_c] = size(S);
n_R = sum(abs(S),2);


alpha = MX.sym('alpha',n_c);
theta = MX.sym('theta',n_sys);
%%
s = 0;
ii = 1;
chi = (1-s)/2+s*alpha(ii);

%%


g_lift_theta = [];
g_lift_beta = [];
beta = [];
% for ii = 1:n_c
        current_beta = 1;
        current_beta_one = 1;           
        current_beta_minus_one = 1;           
        for n_current = 1:n_c
            ind_current = find(n_R >= n_current);
            temp_S = S(ind_current,1:n_current);
            % index w.r.t. temp matrix;
            ind_one = find(temp_S(:,n_current )  == 1);
            ind_minus_one = find(temp_S(:,n_current)  == -1);
            % index w.r.t. temp S matrix;
            ind_one = ind_current(ind_one);
            ind_minus_one = ind_current(ind_minus_one);
            if length(ind_one) == 1
                     eval(['g_lift_theta = [g_lift_theta; theta(ind_one) - current_beta *alpha(ind_one)];'])
                % make theta of this row
            else
            if 1<n_current && n_current<n_c
                % update beta current
                % \beta_{number_of_dive,number_of_lift_variable_in_this_dive}$
                % define variable
                 eval(['beta_' num2str(n_current-1) '1= MX.sym(''beta_' num2str(n_current-1) '1'',1);'])
                 % add to other beta
                 eval(['beta = [beta;beta_' num2str(n_current-1) '1];'])
                 % add to the lifting constraint
                 eval(['g_lift_beta  = [g_lift_beta ; beta_' num2str(n_current-1) '1 - current_beta*alpha(n_current)];'])
                 % update current working beta
                 eval(['current_beta  = beta_' num2str(n_current-1) '1;'])
            else
                current_beta = current_beta*alpha(n_current); 
            end
            
            end

            if length(ind_minus_one) == 1
                % make theta of this row
                     eval(['g_lift_theta = [g_lift_theta; theta(ind_minus_one) - current_beta*(1-alpha(n_current))];'])
            else
                 if 1<n_current && n_current<n_c
                   % update beta current
                % \beta_{number_of_dive,number_of_lift_variable_in_this_dive}$
                % define variable
                 eval(['beta_' num2str(n_current-1) '2= MX.sym(''beta_' num2str(n_current-1) '2'',1);'])
                 % add to other beta
                 eval(['beta = [beta;beta_' num2str(n_current-1) '2];'])
                 % add to the lifting constraint
                 eval(['g_lift_beta  = [g_lift_beta ; beta_' num2str(n_current-1) '2 - current_beta*(1-alpha(n_current))];'])
                 % update current working beta
                 eval(['current_beta = beta_' num2str(n_current-1) '2;'])
                 else
                     current_beta = current_beta*(1-alpha(n_current));
                 end
            end                                   
        end
% end
beta
g_lift_beta
g_lift_theta
if length(g_lift_theta) == n_sys
    disp('sucess;')
else
    disp('fail;')
end