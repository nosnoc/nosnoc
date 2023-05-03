function [J,g_comp,g_comp_lb,g_comp_ub] = reformulate_mpcc_constraints(varargin)
% treatment of the bilinear constraint that arise from complementarity or orthogonality conditions
%   [J,g_comp,g_comp_lb,g_comp_ub] = reformulate_mpcc_constraints(objective_scaling_direct,mpcc_mode,mpcc_var_this_fe,dims,current_index);
% Sanity check of the input
% List of mpcc_modes
% 'direct'
% 'Scholtes_eq'
% 'Scholtes_ineq'
% 'ell_1_penalty'
% 'elastic_ineq'
% 'elastic_eq'
% 'elastic_two_sided'
% 'elastic_ell_1_ineq'
% 'elastic_ell_1_eq'
% 'elastic_ell_1_two_sided'
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
dims = varargin{4};
unfold_struct(dims,'caller')

% TODO: no vargin here!

% Indices
if nargin == 4
    mpcc_mode = 'ell1_penalty';
    k = N_stages-1;
    i = N_finite_elements(end)-1;
    j = d;
else
    current_index = varargin{5};
    %     k = current_index.k;  i = current_index.i; j = current_index.j;
    unfold_struct(current_index,'caller')
end
%% Initialize
n_all_comp_j = length(g_all_comp_j); % outputdimension of the bilinear constraint.
g_comp = [];
g_comp_lb = [];
g_comp_ub = [];
%% Adding all Complementarity Constraints (standard and cross complementarity), their treatment depends on the chosen MPCC Method.
if ~isempty(g_all_comp_j)
    if strcmpi(mpcc_mode,'direct')
        % Info: this is implemented as a special case of scholtes_ineq, with sigma_k = 0';
        % works better with X1x2 <=0 than X1x2 = 0;
    elseif strcmpi(mpcc_mode,'scholtes_eq')
        % smoothed  - algebraic constraints
        g_comp = g_all_comp_j-p(1)*ones(n_all_comp_j,1);
        g_comp_lb = zeros(n_all_comp_j,1);
        g_comp_ub = zeros(n_all_comp_j,1);
    elseif strcmpi(mpcc_mode,'scholtes_ineq')
        g_comp = g_all_comp_j-p(1)*ones(n_all_comp_j,1);
        g_comp_lb = -inf*ones(n_all_comp_j,1);
        g_comp_ub = zeros(n_all_comp_j,1);
    elseif strcmpi(mpcc_mode,'ell_1_penalty')
        if objective_scaling_direct
            J = J + (1/p(1))*sum(g_all_comp_j);
        else
            J = p(1)*J + sum(g_all_comp_j);
        end
    elseif strcmpi(mpcc_mode,'elastic_ineq')
        % elastic mode - inequality
        g_comp = g_all_comp_j-s_elastic*ones(n_all_comp_j,1);
        g_comp_lb = -inf*ones(n_all_comp_j,1);
        g_comp_ub = zeros(n_all_comp_j,1);
    elseif strcmpi(mpcc_mode,'elastic_eq')
        % elastic mode - equality
        g_comp = g_all_comp_j-s_elastic*ones(n_all_comp_j,1);
        g_comp_lb = zeros(n_all_comp_j,1);
        g_comp_ub = zeros(n_all_comp_j,1);
    elseif strcmpi(mpcc_mode,'elastic_two_sided')
        % elastic mode - two sided inequality
        g_comp = g_all_comp_j-s_elastic*ones(n_all_comp_j,1);
        g_comp_ub = zeros(n_all_comp_j,1);
        g_comp_lb = -inf*ones(n_all_comp_j,1);

        g_comp = [g_comp;g_all_comp_j+s_elastic*ones(n_all_comp_j,1)];
        g_comp_ub = [g_comp_ub; inf*ones(n_all_comp_j,1)];
        g_comp_lb = [g_comp_lb;  zeros(n_all_comp_j,1)];
    elseif strcmpi(mpcc_mode,'elastic_ell_1_ineq')
        g_comp = g_all_comp_j-s_elastic;
        g_comp_lb = -inf*ones(n_all_comp_j,1);
        g_comp_ub = zeros(n_all_comp_j,1);
    elseif strcmpi(mpcc_mode,'elastic_ell_1_eq')
        g_comp = g_all_comp_j-s_elastic;
        g_comp_lb = zeros(n_all_comp_j,1);
        g_comp_ub = zeros(n_all_comp_j,1);
    elseif strcmpi(mpcc_mode,'elastic_ell_1_two_sided')
        % elastic mode - two sided inequality
        g_comp = g_all_comp_j-s_elastic;
        g_comp_ub = zeros(n_all_comp_j,1);
        g_comp_lb = -inf*ones(n_all_comp_j,1);

        g_comp = [g_comp;g_all_comp_j+s_elastic];
        g_comp_ub = [g_comp_ub; inf*ones(n_all_comp_j,1)];
        g_comp_lb = [g_comp_lb;  zeros(n_all_comp_j,1)];
    else
        error('Pick a valid option for mpcc_mode, e.g., ''Scholtes_ineq''.');
    end
end
end

