function [J,g_comp,g_comp_lb,g_comp_ub] = reformulate_mpcc_constraints(varargin)
% treatment of the bilinear constraint that arise from complementarity or orthogonality conditions
%   [J,g_comp,g_comp_lb,g_comp_ub] = reformulate_mpcc_constraints(objective_scaling_direct,mpcc_mode,mpcc_var_this_fe,dimensions,current_index);


% Sanity check of the input
import casadi.*
%%
% g_cross_comp_j = [];
% Sanity Chack
if ~(nargin == 4 || nargin == 5)
    error('Wrong input number, see function help for required inputs.')
end

% Assign variable names to all inputs
% FESD settings
objective_scaling_direct = varargin{1};
mpcc_mode = varargin{2};
% Complementarity variables
mpcc_var_this_fe = varargin{3};
unfold_struct(mpcc_var_this_fe,'caller')

% Dimensions
dimensions = varargin{4};
unfold_struct(dimensions,'caller')

% Indices
if nargin == 4
    mpcc_mode = 4;
    k = N_stages-1;
    i = N_finite_elements(end)-1;
    j = d;
else
    current_index = varargin{5};
%     k = current_index.k;  i = current_index.i; j = current_index.j;
    unfold_struct(current_index,'caller')
end
%% Initalize
n_all_comp_j = length(g_all_comp_j); % outputdimension of the bilinear constraint.
 g_comp = [];
 g_comp_lb = [];
 g_comp_ub = [];
 %% Adding all Complementarity Constraints (standard and cross complementarity), their treatment depends on the chosen MPCC Method.
             if ~isempty(g_all_comp_j)
                switch mpcc_mode
                    case 1
                        disp('Info: this is implemented as a special case of 3, with sigma = 0');
                        % works better with X1x2 <=0 than X1x2 = 0;
                    case 2
                        % smoothed  - algebraic constraints
                        g_comp = g_all_comp_j-p*ones(n_all_comp_j,1);
                        g_comp_lb = zeros(n_all_comp_j,1);
                        g_comp_ub = zeros(n_all_comp_j,1);
                    case 3
                        g_comp = g_all_comp_j-p*ones(n_all_comp_j,1);
                        g_comp_lb = -inf*ones(n_all_comp_j,1);
                        g_comp_ub = zeros(n_all_comp_j,1);
                    case 4
                        if k == N_stages-1 && j == n_s && i == N_finite_elements(end)-1
%                             if cross_comp_mode == 10 &&  k == N_stages-1 && j == d && i == N_finite_elements(end)-1
                            if objective_scaling_direct
                                J = J + (1/p)*g_all_comp_j;
                            else
                                J = p*J + g_all_comp_j;
                            end
                        end
                    case 5
                        % elastic mode - inequality
                        g_comp = g_all_comp_j-s_elastic*ones(n_all_comp_j,1);
                        g_comp_lb = -inf*ones(n_all_comp_j,1);
                        g_comp_ub = zeros(n_all_comp_j,1);
                    case 6
                        % elastic mode - equality
                        g_comp = g_all_comp_j-s_elastic*ones(n_all_comp_j,1);
                        g_comp_lb = zeros(n_all_comp_j,1);
                        g_comp_ub = zeros(n_all_comp_j,1);
                    case 7
                        % elastic mode - two sided inequality
                        g_comp = g_all_comp_j-s_elastic*ones(n_all_comp_j,1);
                        g_comp_ub = zeros(n_all_comp_j,1);
                        g_comp_lb = -inf*ones(n_all_comp_j,1);

                        g_comp = [g_comp;g_all_comp_j+s_elastic*ones(n_all_comp_j,1)];
                        g_comp_ub = [g_comp_ub; inf*ones(n_all_comp_j,1)];
                        g_comp_lb = [g_comp_lb;  zeros(n_all_comp_j,1)];
                        
                    case 8
                        % elastic mode with barrier homotopy
                        % elastic mode - inequality
                        g_comp = g_all_comp_j-s_elastic*ones(n_all_comp_j,1);
                        g_comp_lb = -inf*ones(n_all_comp_j,1);
                        g_comp_ub = zeros(n_all_comp_j,1);
%                         lbg = [lbg; g_comp_lb];
%                         ubg = [ubg; g_comp_ub];
                    case 9
                        % elastic mode with barrier homotopy
                        % elastic mode - equality
                        g_comp = g_all_comp_j-s_elastic*ones(n_all_comp_j,1);
                        g_comp_lb = zeros(n_all_comp_j,1);
                        g_comp_ub = zeros(n_all_comp_j,1);

                    case 10
                        if  k == N_stages-1 && j == n_s && i == N_finite_elements(end)-1
                            if objective_scaling_direct
                                J = J + (1/s_elastic)*g_all_comp_j;
                            else
                                J = s_elastic*J + g_all_comp_j;
                            end
                        end
                    case 11
                          % ell_1 elastic mode - inequality
                        g_comp = g_all_comp_j-s_elastic;
                        g_comp_lb = -inf*ones(n_all_comp_j,1);
                        g_comp_ub = zeros(n_all_comp_j,1);
                        case 12
                       % elastic mode - equality
                        g_comp = g_all_comp_j-s_elastic;
                        g_comp_lb = zeros(n_all_comp_j,1);
                        g_comp_ub = zeros(n_all_comp_j,1);
                    case 13
                        % elastic mode - two sided inequality
                        g_comp = g_all_comp_j-s_elastic;
                        g_comp_ub = zeros(n_all_comp_j,1);
                        g_comp_lb = -inf*ones(n_all_comp_j,1);

                        g_comp = [g_comp;g_all_comp_j+s_elastic];
                        g_comp_ub = [g_comp_ub; inf*ones(n_all_comp_j,1)];
                        g_comp_lb = [g_comp_lb;  zeros(n_all_comp_j,1)];
                        

                    otherwise
                        error('Pick mpcc_mode between 1 and 13.');
                end
            end

end

