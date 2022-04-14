function [J_comp,g_cross_comp_j] = create_complementarity_constraints_old(varargin)
% A function for formulating the complementarity and cross % complementarity(orthogonaliy constaints)
% Examples of calling this function
% [J_comp,g_cross_comp_j] = create_complementarity_constraints(use_fesd,cross_comp_mode,comp_var_current_fe,dimensions,current_index);
% [J_comp,g_cross_comp_j] = create_complementarity_constraints(use_fesd,cross_comp_mode,comp_var_current_fe,dimensions);

import casadi.*
%%
g_cross_comp_j = [];
g_cross_comp_temp = 0;
% Sanity Chack
if ~(nargin == 4 || nargin == 5)
    error('Wrong input number, see function help for required inputs.')
end

% Assign variable names to all inputs
% FESD settings
use_fesd = varargin{1};
cross_comp_mode = varargin{2};
% Complementarity variables
comp_var_current_fe = varargin{3};
unfold_struct(comp_var_current_fe,'caller')

% Dimensions
dimensions = varargin{4};
unfold_struct(dimensions,'caller')

% Indices
if nargin == 4
    cross_comp_mode = 10;
    k = N_stages-1;
    i = N_finite_elements(end)-1;
    j = n_s;
else
    current_index = varargin{5};
    unfold_struct(current_index,'caller')
end
%% Formlation of standard and cross complementarity consraints.
if use_fesd
    % Formulae for the cross complementarities
    g_cross_comp_j = [];
    for r = 1:n_simplex
        % for different subsystems the lambdas and thetas are decoupled and should be treated as such.
        ind_temp_theta = m_ind_vec(r):m_ind_vec(r)+m_vec(r)-1;
        ind_temp_lambda = ind_temp_theta;

        Theta_ki_j = Theta_ki_current_fe{j}(ind_temp_theta);

        % sum of all cross-complementarities (vector-valued) --> later put into scalar value for the complementarity residual
        J_comp  = J_comp + diag(Theta_ki_j)*sum_lambda_ki(ind_temp_theta);

        switch cross_comp_mode
            case 1
                % Full Sparsity. Every point with every gives a vector valued constraint
                if i > 0
                    % if it is not the very first element, include lambda_{n,q,0} (\lambda_{n,0})
                    Lambda_ki_jj = Lambda_end_previous_fe{1}(ind_temp_lambda);
                    g_cross_comp_j = [g_cross_comp_j ;diag(Theta_ki_j)*(Lambda_ki_jj)];
                end
                for jj = 1:n_s
                    Lambda_ki_jj = Lambda_ki_current_fe{jj}(ind_temp_lambda);
                    g_cross_comp_j = [g_cross_comp_j ;diag(Theta_ki_j)*(Lambda_ki_jj)];
                end
            case 2
                % Full Sparsity. Every point with every gives a scalar valued constraint
                if i > 0
                    % if it is not the very first element, include lambda_{n,q,0} (\lambda_{n,0})
                    Lambda_ki_jj = Lambda_end_previous_fe{1}(ind_temp_lambda);
                    g_cross_comp_j = [g_cross_comp_j ;(Theta_ki_j)'*(Lambda_ki_jj)];
                end
                for jj = 1:n_s
                    Lambda_ki_jj = Lambda_ki_current_fe{jj}(ind_temp_lambda);
                    g_cross_comp_j = [g_cross_comp_j ;(Theta_ki_j)'*(Lambda_ki_jj)];
                end

                % Casese 3 and 4: for every stage point a vector valued constraint
            case 3
                %  For every stage point one vector-valued  constraint via sum of all \lambda
                g_cross_comp_j = [g_cross_comp_j ;diag(Theta_ki_j)*sum_lambda_ki(ind_temp_theta)];
            case 4
                %  Case 4: For every stage point one vector-valued constraint via sum of all \theta.
                if j == 1 && i > 0
                    % take care of lambda_{n,0} at same time as \lambda_{n,1}
                    Lambda_ki_jj = Lambda_end_previous_fe{1}(ind_temp_lambda);
                    g_cross_comp_j = [g_cross_comp_j;diag(Lambda_ki_jj)*sum_theta_ki];
                end
                Lambda_ki_jj = Lambda_ki_current_fe{j}(ind_temp_lambda);
                g_cross_comp_j = [g_cross_comp_j;diag(Lambda_ki_jj)*sum_theta_ki(ind_temp_theta)];

                % Cases 5 and 6: For every stage point a scalarv alued constraint
            case 5
                % Case 5: Per stage point one scalar constraint via sum of \lambda.
                g_cross_comp_j = [g_cross_comp_j ;(Theta_ki_j)'*sum_lambda_ki(ind_temp_theta)];
            case 6
                %  Case 6: For every stage point one scalar-valued constraint via sum of all \theta.
                if j == 1 && i > 0
                    % take care of lambda_{n,0} at same time as \lambda_{n,1}
                    Lambda_ki_jj = Lambda_end_previous_fe{1}(ind_temp_lambda);
                    g_cross_comp_j = [g_cross_comp_j;(Lambda_ki_jj)'*sum_theta_ki];
                end
                Lambda_ki_jj = Lambda_ki_current_fe{j}(ind_temp_lambda);
                g_cross_comp_j = [g_cross_comp_j;(Lambda_ki_jj)'*sum_theta_ki(ind_temp_theta)];

                % Cases 7 and 8: a vector or single valued constraint per every finite element.
            case 7
                %  constraint via sum of \lambda
                g_cross_comp_temp = g_cross_comp_temp +diag(Theta_ki_j)*sum_lambda_ki(ind_temp_theta);
                if j == n_s
                    g_cross_comp_j = [g_cross_comp_j ;g_cross_comp_temp];
                end
            case 8
                %  Case 8: Per stage (finite element) one scalar constraint via sum of \lambda.
                if j == n_s
                    % add it only once during the loop
                    g_cross_comp_j = [g_cross_comp_j ;(sum_lambda_ki)'*(sum_theta_ki)];
                end

                % cases 9 and 10 are fully integral options either vector or scalar valued.
            case 9
                if k == N_stages-1 && j == n_s && i == N_finite_elements(end)-1
                    g_cross_comp_j = [g_cross_comp_j ;J_comp];
                end
            case 10
                if k == N_stages-1 && j == n_s && i == N_finite_elements(end)-1
                    g_cross_comp_j = [g_cross_comp_j ;sum(J_comp)];
                end
            otherwise
                error('Please pick cross_comp_mode between 1 and 10.')
        end
    end
else
    g_cross_comp_j = diag(Lambda_ki_current_fe{j})*Theta_ki_current_fe{j};
    J_comp = J_comp+sum(g_cross_comp_j);
end
end

