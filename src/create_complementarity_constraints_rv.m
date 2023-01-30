% BSD 2-Clause License

% Copyright (c) 2022, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% This file is part of NOSNOC.

%
%
function [results_cross_comp] = create_complementarity_constraints_rv(varargin)

% A function for formulating the complementarity and cross % complementarity(orthogonaliy constaints)
% Examples of calling this function
% [J_comp,g_cross_comp_j] = create_complementarity_constraints(use_fesd,cross_comp_mode,comp_var_current_fe,dimensions,current_index);
% [J_comp,g_cross_comp_j] = create_complementarity_constraints(use_fesd,cross_comp_mode,comp_var_current_fe,dimensions);
import casadi.*
%% Process input
% basic
use_fesd = varargin{1};
cross_comp_mode = varargin{2};
comp_var_current_fe = varargin{3};
unfold_struct(comp_var_current_fe,'caller')
% Dimensions
dimensions = varargin{4};
unfold_struct(dimensions,'caller')
% MPCC Function & parameters
Psi_mpcc_fun = varargin{5};
sigma_p = varargin{6};
% Indices
if nargin == 6
    cross_comp_mode = 10;
    k = N_stages-1;
    i = N_finite_elements(end)-1;
    j = n_s;
else
    current_index = varargin{7};
    unfold_struct(current_index,'caller');

end
% Empty vectors for complementarity constraints
g_cross_comp_ki = [];

%% Formlation of standard and cross complementarity consraints.
n_lambda = size(Lambda_ki,2);% ; number of lambda_ki in finite element
% main loop;
for r = 1:n_sys
    % for different subsystems the lambdas and thetas are decoupled and should be treated as such.
    ind_temp_theta = m_ind_vec(r):m_ind_vec(r)+m_vec(r)-1;
    Theta_ki_temp = Theta_ki(ind_temp_theta,:);
    Lambda_ki_temp = Lambda_ki(ind_temp_theta,:);
    if use_fesd

        %% update vector valued sumes over control interval, update only once per finite element
        g_comp = zeros(n_theta,1);
        g_comp_residual = 0;
        % sum for current finite element
        for ii = 1:n_s
            for jj = 1:n_lambda
                if ~is_zero(Lambda_ki_temp(:,jj))
                    g_comp  = g_comp + Psi_mpcc_fun(Theta_ki_temp(:,ii),Lambda_ki_temp(:,jj),sigma_p);
                    g_comp_residual = g_comp_residual+Theta_ki_temp(:,ii)'*Lambda_ki_temp(:,jj);
                end
            end
        end
        % MPCC function updates (for constraints)
        g_cross_comp_k = g_cross_comp_k + g_comp;
        g_cross_comp_all = g_cross_comp_all + g_comp;
        % Pure bilinear updates (for complementarity evaluation)
        cross_comp_residual_k = cross_comp_residual_k + g_comp_residual;
        cross_comp_residual_all = cross_comp_residual_all  + g_comp_residual;



        %% Complementarity constraint formulation
        switch cross_comp_mode
            % Cases 1 and 2: Full sparsity,
            case 1
                %                  n_theta*n_s*(n_s+1)  constraints pre finite element
                for jj = 1:n_lambda
                    g_comp  = Psi_mpcc_fun(Theta_ki_temp(:),repmat((Lambda_ki_temp(:,jj)),n_s,1),sigma_p);
                    g_cross_comp_ki = [g_cross_comp_ki;g_comp];
                end
            case 2
                % Full Sparsity. Every point with every gives a scalar-valued constraint, n_s*(n_s+1)  constraints pre finite element
                for ii = 1:n_s
                    for jj = 1:n_lambda
                        if is_zero(Lambda_ki_temp(:,jj))
                            g_comp   = [];
                        else
                            g_comp  = Psi_mpcc_fun(Theta_ki_temp(:,ii),Lambda_ki_temp(:,jj),sigma_p/(n_s+1));
                        end
                        g_cross_comp_ki = [g_cross_comp_ki; sum(g_comp)];
                    end
                end
                % Cases 3 to 6: For every stage j one constraint, 3 and 4 vector valued, 5 and 6 scalar valued
            case 3
                for ii = 1:n_s
                    g_comp = zeros(n_theta,1);
                    for jj = 1:n_lambda
                        g_comp  = g_comp + Psi_mpcc_fun(Theta_ki_temp(:,ii),(Lambda_ki_temp(:,jj)),sigma_p/(n_s+1));
                    end
                    g_cross_comp_ki = [g_cross_comp_ki;g_comp];
                end
            case 4
                %  Case 4: For every stage point one vector-valued constraint via sum of all \theta.
                for jj = 1:n_lambda
                    g_comp = zeros(n_theta,1);
                    for ii = 1:n_s
                        if is_zero(Lambda_ki_temp(:,jj))
                            g_comp   = [];
                        else
                            g_comp  = g_comp + Psi_mpcc_fun(Theta_ki_temp(:,ii),Lambda_ki_temp(:,jj),sigma_p/(n_s));
                        end
                    end
                    g_cross_comp_ki = [g_cross_comp_ki;g_comp];
                end
            case 5
                % Case 5: Per stage point one scalar constraint via sum of  \lambda.
                % n_s constraints, same as 3 but with summing all terms
                for ii = 1:n_s
                    g_comp = zeros(n_theta,1);
                    for jj = 1:n_lambda
                        g_comp  = g_comp + Psi_mpcc_fun(Theta_ki_temp(:,ii),(Lambda_ki_temp(:,jj)),sigma_p/(n_s));
                    end
                    g_cross_comp_ki = [g_cross_comp_ki;sum(g_comp)];
                end
            case 6
                %  Case 6: For every stage point one scalar-valued constraint via sum of all \theta.
                for jj = 1:n_lambda
                    g_comp = zeros(n_theta,1);
                    for ii = 1:n_s
                        if is_zero(Lambda_ki_temp(:,jj))
                            g_comp   = [];
                        else
                            g_comp  = g_comp + Psi_mpcc_fun(Theta_ki_temp(:,ii),Lambda_ki_temp(:,jj),sigma_p/(n_s));
                        end
                    end
                    g_cross_comp_ki = [g_cross_comp_ki;sum(g_comp)];
                end
                % Cases 7 and 8: a vector or scalar valued constraint per every finite element
            case 7
                % n_theta constraint per finite elements
                g_comp = zeros(n_theta,1);
                for ii = 1:n_s
                    for jj = 1:n_lambda
                        %                         g_comp  = g_comp + Psi_mpcc_fun(Theta_ki_temp(:,ii),(Lambda_ki_temp(:,jj)),sigma_p/(n_s*(n_s+1)));
                        g_comp  = g_comp + Psi_mpcc_fun(Theta_ki_temp(:,ii),(Lambda_ki_temp(:,jj)),sigma_p);
                    end
                end
                g_cross_comp_ki = [g_cross_comp_ki;g_comp];
            case 8
                % same as 7, but scalar valued, one scalar valed per finite elemetn
                g_comp = zeros(n_theta,1);
                for ii = 1:n_s
                    for jj = 1:n_lambda
                        %                         g_comp  = g_comp + Psi_mpcc_fun(Theta_ki_temp(:,ii),(Lambda_ki_temp(:,jj)),sigma_p/(n_s*(n_s+1)));
                        g_comp  = g_comp + Psi_mpcc_fun(Theta_ki_temp(:,ii),(Lambda_ki_temp(:,jj)),sigma_p);
                    end
                end
                g_cross_comp_ki = [g_cross_comp_ki;sum(g_comp)];
                %% Cases 9 and 10, per control intevral one constraint
            case 9
                % one vector-valued per control interval (n_theta constraints)
                if i == N_finite_elements(k+1)-1
                    g_cross_comp_ki = [g_cross_comp_ki ;g_cross_comp_k];
                end
            case 10
                % one scalar valued per control interval
                if i == N_finite_elements(k+1)-1
                    % the if is because the sum over control interval makes only sense at this point
                    g_cross_comp_ki = [g_cross_comp_ki ;sum(g_cross_comp_k)];
                end
                %  case 11 and 12: single constraint over whole horizon
            case 11
                % one vector-valued for all control intervals
                if k == N_stages-1 && i == N_finite_elements(end)-1
                    g_cross_comp_ki = [g_cross_comp_ki ;g_cross_comp_all];
                end
            case 12
                % one scalar vector-valued for all control intervals
                if k == N_stages-1 && i == N_finite_elements(end)-1
                    g_cross_comp_ki = [g_cross_comp_ki ;sum(g_cross_comp_all)];
                end
            otherwise
                error('Please pick cross_comp_mode between 1 and 12.')
        end
    else
        % TODO: add some sparsity control
        switch cross_comp_mode
            case 1
                g_comp  = Psi_mpcc_fun(Theta_ki_temp(:),Lambda_ki_temp(:),sigma_p);
                g_cross_comp_ki = [g_cross_comp_ki;g_comp];
            case 2
                g_comp  = Psi_mpcc_fun(Theta_ki_temp(:),Lambda_ki_temp(:),sigma_p);
                g_cross_comp_ki = [g_cross_comp_ki;g_comp];
            case 3

            otherwise
        end

        %     g_cross_comp_ki = sum(Psi_mpcc_fun(Lambda_ki(:,2:end)*Theta_ki(:),sigma_p));
    end
end
%% Output
% g_cross_comp_ki,n_cross_comp_ki
results_cross_comp.g_cross_comp_ki = g_cross_comp_ki;
results_cross_comp.n_cross_comp_ki = length(g_cross_comp_ki);
% Constraints
results_cross_comp.g_cross_comp_k = g_cross_comp_k;
results_cross_comp.g_cross_comp_all = g_cross_comp_all;
% Cross-complementarity residual
results_cross_comp.cross_comp_residual_k = cross_comp_residual_k;
results_cross_comp.cross_comp_residual_all = cross_comp_residual_all;
end