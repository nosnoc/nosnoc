% Copyright 2022 Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% This file is part of NOSNOC.

% The 2-Clause BSD License

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
%
function [results_cross_comp] = create_complementarity_constraints(varargin)

% A function for formulating the complementarity and cross complementarity (orthogonaliy constaints)
% Examples of calling this function
% [J_comp,g_cross_comp_j] = create_complementarity_constraints(use_fesd,cross_comp_mode,comp_var_current_fe,dimensions,current_index);
% [J_comp,g_cross_comp_j] = create_complementarity_constraints(use_fesd,cross_comp_mode,comp_var_current_fe,dimensions);
import casadi.*
%%
g_cross_comp_j = [];
% Sanity Chack
if ~(nargin == 4 || nargin == 5)
    error('Wrong input number, see function help for required inputs.')
end
% Assign variable names to all inputs
% FESD settings
use_fesd = varargin{1};
cross_comp_mode = varargin{2};
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
    unfold_struct(current_index,'caller');
end

%% Formulation of standard and cross complementarity consraints.
if use_fesd
    % Formulate for the cross complementarities
    g_cross_comp_j = [];
    % update vector valued sumes over control interval
    if j == n_s 
        % update only once per finite element
        cross_comp_k = cross_comp_k + diag(sum_Theta_ki)*sum_Lambda_ki;
        cross_comp_all = cross_comp_all + diag(sum_Theta_ki)*sum_Lambda_ki;
    end
    for r = 1:n_simplex
        % for different subsystems the lambdas and thetas are decoupled and should be treated as such.
        ind_temp_theta = m_ind_vec(r):m_ind_vec(r)+m_vec(r)-1;
        ind_temp_lambda = ind_temp_theta;
        try
            Theta_kij = Theta_ki{j}(ind_temp_theta);
        catch
            Theta_kij = Theta_ki(ind_temp_theta,j);
        end
        % sum of all cross-complementarities (vector-valued) --> later put into scalar value for the complementarity residual
%         J_comp  = J_comp + diag(Theta_kij)*sum_Lambda_ki(ind_temp_theta);

        switch cross_comp_mode
            %% Cases 1 and 2: Full sparsity, ~ n_s*(n_s+1) constraints
            case 1
                % Full Sparsity. Every point with every gives a vector valued constraint
                if i > 0
                    % if it is not the very first element, include lambda_{n,q,0} (\lambda_{n,0})
                    Lambda_kijj = Lambda_end_previous_fe(ind_temp_lambda);
                    g_cross_comp_j = [g_cross_comp_j;diag(Theta_kij)*(Lambda_kijj)];
                end

                for jj = 1:n_s
                    Lambda_kijj = Lambda_ki{jj}(ind_temp_lambda);
                    g_cross_comp_j = [g_cross_comp_j ;diag(Theta_kij)*(Lambda_kijj)];
                end
            case 2
                % Full Sparsity. Every point with every gives a scalar valued constraint
                if i > 0
                    % if it is not the very first element, include lambda_{n,q,0} (\lambda_{n,0})
                    Lambda_kijj = Lambda_end_previous_fe(ind_temp_lambda);
                    g_cross_comp_j = [g_cross_comp_j ;(Theta_kij)'*(Lambda_kijj)];
                end
                for jj = 1:n_s
                    try
                    Lambda_kijj = Lambda_ki{jj}(ind_temp_lambda);
                    catch
                        Lambda_kijj = Lambda_ki(ind_temp_lambda,jj);
                    end
                    g_cross_comp_j = [g_cross_comp_j ;(Theta_kij)'*(Lambda_kijj)];
                end

                %% Cases 3 to 6: For every stage j one constraint, 3 and 4 vector valued, 5 and 6 scalar valued
            case 3
                %  For every stage point one vector-valued  constraint via sum of all \lambda
                g_cross_comp_j = [g_cross_comp_j ;diag(Theta_kij)*sum_Lambda_ki(ind_temp_theta)];
            case 4
                %  Case 4: For every stage point one vector-valued constraint via sum of all \theta.
                if j == 1 && i > 0
                    % take care of lambda_{n,0} at same time as \lambda_{n,1}
                    Lambda_kijj = Lambda_end_previous_fe(ind_temp_lambda);
                    g_cross_comp_j = [g_cross_comp_j;diag(Lambda_kijj)*sum_Theta_ki];
                end
                Lambda_kijj = Lambda_ki{j}(ind_temp_lambda);
                g_cross_comp_j = [g_cross_comp_j;diag(Lambda_kijj)*sum_Theta_ki(ind_temp_theta)];
                % Cases 5 and 6: For every stage point a scalarv alued constraint
            case 5
                % Case 5: Per stage point one scalar constraint via sum of \lambda.
                g_cross_comp_j = [g_cross_comp_j ;(Theta_kij)'*sum_Lambda_ki(ind_temp_theta)];
            case 6
                %  Case 6: For every stage point one scalar-valued constraint via sum of all \theta.
                if j == 1 && i > 0
                    % take care of lambda_{n,0} at same time as \lambda_{n,1}
                    Lambda_kijj = Lambda_end_previous_fe(ind_temp_lambda);
                    g_cross_comp_j = [g_cross_comp_j;sum_Theta_ki(ind_temp_theta)'*Lambda_kijj(ind_temp_theta)];
                end
                Lambda_kijj = Lambda_ki{j}(ind_temp_lambda);
                g_cross_comp_j = [g_cross_comp_j;sum_Theta_ki(ind_temp_theta)'*Lambda_kijj(ind_temp_theta)];

                %% Cases 7 and 8: a vector or scalar valued constraint per every finite element
            case 7
                if j == n_s
                    g_cross_comp_j = [g_cross_comp_j ;diag(sum_Theta_ki(ind_temp_theta))*sum_Lambda_ki(ind_temp_theta)];
                end
            case 8
                %  Case 8: same as 7, but vector valued
                if j == n_s
                    temp = diag(sum_Theta_ki(ind_temp_theta))*sum_Lambda_ki(ind_temp_theta);
                    g_cross_comp_j = [g_cross_comp_j ;sum(temp)];
                end
                %% Cases 9 and 10, per control intevral one constraint
            case 9
                % vectorvalued
                if i == N_finite_elements(k+1)-1 && j == n_s
                    % the if is because the sum over control interval makes only sense at this point
                    g_cross_comp_j = [g_cross_comp_j ;cross_comp_k];
                end
            case 10
                % vectorvalued
                if i == N_finite_elements(k+1)-1 && j == n_s
                    % the if is because the sum over control interval makes only sense at this point
                    g_cross_comp_j = [g_cross_comp_j ;sum(cross_comp_k)];
                end
                %%  case 11 and 12: ssingle constraint over whole horizon
             case 11
                 % vector valued
                if k == N_stages-1 && j == n_s && i == N_finite_elements(end)-1
                    g_cross_comp_j = [g_cross_comp_j ;cross_comp_all];
                end
            case 12
                % scalar valued
                if k == N_stages-1 && j == n_s && i == N_finite_elements(end)-1
                    g_cross_comp_j = [g_cross_comp_j ;sum(cross_comp_all)];
                end
            otherwise
                error('Please pick cross_comp_mode between 1 and 12.')
        end
    end
else
    g_cross_comp_j = diag(Lambda_ki{j})*Theta_ki{j};
end
%% Output
results_cross_comp.g_cross_comp_j = g_cross_comp_j;
results_cross_comp.cross_comp_k = cross_comp_k;
results_cross_comp.cross_comp_all = cross_comp_all;

end

