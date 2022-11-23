%    This file is part of NOSNOC.
%
%    NOSNOC -- A software for NOnSmooth Numerical Optimal Control.
%    Copyright (C) 2022 Armin Nurkanovic, Moritz Diehl (ALU Freiburg).
%
%    NOSNOC is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 3 of the License, or (at your option) any later version.
%
%    NOSNOC is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with NOSNOC; if not, write to the Free Software Foundation,
%    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%
% function [solver,solver_initialization, model,settings] = create_nlp_nosnoc(model,settings)
function [varargout] = create_nlp_nosnoc(varargin)
% This functions creates the solver instance for the OCP discretized with FESD (or time-stepping IRK scheme).
% The discretization results in an MPCC which can be solved by various
% reformulations, see below.
% -------------------
% Brief MPCC Wiki
% There are several possible MPCC Solution strategies avilable, by setting mpcc_mode to :
% 1 - treat complementarity conditions directly in the NLP, the bilinear term is tread as an inequality constraint.
% 2 - Smooth the complementarity conditions.
% 3 - Relax the complementarity conditions.
% 4 - \ell_1 penalty, penalize the sum of all bilinear terms in the objective
% 5 - \ell_infty elastic mode, upper bound all bilinear term with a positive slack, and penalize the slack in the objective.
% 6 - \ell_infty elastic mode, equate all bilinear term to a positive slack, and penalize the slack in the objective.
% 7 - \ell_infty, same as 5 but two sided.
% 8 - \ell_infty, same as 5 but the penalty/slack is controlled via the barier formulation
% 9 - \ell_infty, same as 6 but the penalty/slack is controlled via the barier formulation
% 10 - \ell_1, same as 4 but the penalty/slack is controlled via the barier formulation
%% Import CasADi in the workspace of this function
import casadi.*
%% Read data
model = varargin{1};
settings = varargin{2};
%% Reformulation of the PSS into a DCS
[settings] = refine_user_settings(settings);
[model,settings] = model_reformulation_nosnoc(model,settings);

%% Fillin missing settings with default settings
[settings] = fill_in_missing_settings(settings,model);

%% Load user settings and model details
unfold_struct(settings,'caller')
unfold_struct(model,'caller');

%% Parameters
p_val = [sigma_0,rho_sot,rho_h,rho_terminal,T];
% define parameters;
sigma_p = define_casadi_symbolic(casadi_symbolic_mode,'sigma_p'); % homotopy parameter
rho_sot_p = define_casadi_symbolic(casadi_symbolic_mode,'rho_sot_p'); % homotopy parameter
rho_h_p = define_casadi_symbolic(casadi_symbolic_mode,'rho_h_p'); % homotopy parameter
rho_terminal_p = define_casadi_symbolic(casadi_symbolic_mode,'rho_terminal_p'); % homotopy parameter
T_ctrl_p  = define_casadi_symbolic(casadi_symbolic_mode,'T_ctrl_p'); % homotopy parameter

lambda_00 = define_casadi_symbolic(casadi_symbolic_mode,'lambda_00', n_theta);

p = vertcat(sigma_p,rho_sot_p,rho_h_p,rho_terminal_p,T_ctrl_p, lambda_00);

%% Initialization and bounds for step-size
if use_fesd
    ubh = (1+gamma_h)*h_k;
    lbh = (1-gamma_h)*h_k;
    if time_rescaling && ~use_speed_of_time_variables
        % if only time_rescaling is true, speed of time and step size all lumped together, e.g., \hat{h}_{k,i} = s_n * h_{k,i}, hence the bounds need to be extended.
        ubh = (1+gamma_h)*h_k*s_sot_max;
        lbh = (1-gamma_h)*h_k/s_sot_min;
    end
    % initial guess for the step-size
    h0_k = h_k.*ones(N_stages,1);
end

%%  Butcher Tableu
switch irk_representation
    case 'integral'
        [B,C,D,tau_root] = generate_butcher_tableu_integral(n_s,irk_scheme);
        if tau_root(end) == 1
            right_boundary_point_explicit  = 1;
        else
            right_boundary_point_explicit  = 0;
        end
    case 'differential'
        [A_irk,b_irk,c_irk,order_irk] = generate_butcher_tableu(n_s,irk_scheme);
        if c_irk(end) <= 1+1e-9 && c_irk(end) >= 1-1e-9
            right_boundary_point_explicit  = 1;
        else
            right_boundary_point_explicit  = 0;
        end
    otherwise
        error('Choose irk_representation either: ''integral'' or ''differential''')
end

%% Time optimal control
if time_optimal_problem
    % the final time in time optimal control problems
    T_final = define_casadi_symbolic(casadi_symbolic_mode,'T_final',1);
    T_final_guess = T;
end

%% Elastic Mode Variables
s_ell_inf_elastic_exists  = 0;
if mpcc_mode >= 5 && mpcc_mode <= 10
    s_elastic = define_casadi_symbolic(casadi_symbolic_mode,'s_elastic',1);
    s_ell_inf_elastic_exists = 1;
end
s_ell_1_elastic_exists  = 0;

sum_s_elastic = 0;
if mpcc_mode >= 11 && mpcc_mode <= 13
    s_ell_1_elastic_exists  = 1;
end

%% Formulate NLP - Start with an empty NLP
% degrees of freedom
w = {};
w0 = [];
lbw = [];
ubw = [];
% objective
J = 0;
J_comp = 0;
J_comp_std = 0;
J_comp_fesd = 0;
J_regularize_h = 0;
J_regularize_sot = 0;
% constraints
g = {};
lbg = [];
ubg = [];

X_ki = define_casadi_symbolic(casadi_symbolic_mode,'X0',n_x);
w = {w{:}, X_ki};

if there_exist_free_x0
    x0_ub = x0;
    x0_lb = x0;
    x0_ub(ind_free_x0) = inf;
    x0_lb(ind_free_x0) = -inf;

    x0_ub(ind_free_x0) = ubx(ind_free_x0);
    x0_lb(ind_free_x0) = lbx(ind_free_x0);
    lbw = [lbw; x0_lb];
    ubw = [ubw; x0_ub];
else
    lbw = [lbw; x0];
    ubw = [ubw; x0];
end

w0 = [w0; x0];

% Index vectors
ind_x = [1:n_x];
ind_u = [];
ind_v = [];
ind_z = [];
ind_h = [];
ind_elastic = [];
ind_g_clock_state = [];
ind_vf  = [];
ind_sot = []; % index for speed of time variable
ind_boundary = []; % index of bundary value lambda and mu
ind_total = ind_x;

% Integral of the clock state if no time rescaling is present.
sum_h_ki_control_interval_k = 0;
sum_h_ki_all = 0;
% Initialization of forward and backward sums for the step-equilibration
sigma_theta_B_k = 0; % backward sum of theta at finite element k
sigma_lambda_B_k = 0; % backward sum of lambda at finite element k
sigma_lambda_F_k = 0; % forward sum of lambda at finite element k
sigma_theta_F_k = 0; % forward sum of theta at finite element k
nu_vector = [];
% Integral of the clock state if no time-freezing is present.
integral_clock_state = 0;
integral_clock_state_k = 0;
n_cross_comp = zeros(max(N_finite_elements),N_stages);

% initialization
Lambda_end_previous_fe = lambda_00;
Z_kd_end = zeros(n_z,1);
% initialize cross comp and mpcc related structs
mpcc_var_current_fe.p = p;
comp_var_current_fe.cross_comp_k = 0;
comp_var_current_fe.cross_comp_all = 0;

%% Formulate the NLP / Main Discretization loop
for k=0:N_stages-1
    %% Variables for the controls
    if n_u > 0
        Uk = define_casadi_symbolic(casadi_symbolic_mode,['U_' num2str(k)],n_u);
        w = {w{:}, Uk};
        % intialize contros, lower and upper bounds
        w0 = [w0; u0];
        lbw = [lbw;lbu];
        ubw = [ubw;ubu];
        % index colector for contorl variables
        ind_u = [ind_u,ind_total(end)+1:ind_total(end)+n_u];
        ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_u];
    end

    %%  Time rescaling of the stages (speed of time) to acchieve e.g., a desired final time in Time-Freezing or to solve time optimal control problems.
    %     If only time_rescaling is true, then the h_k also make sure to addapt the length of the subintervals, if both
    %     time_rescaling && use_speed_of_time_variables are true, new control variables are introduecd, they can be per stage or one for the whole interval.
    if time_rescaling && use_speed_of_time_variables
        if local_speed_of_time_variable
            % at every stage
            s_sot_k = define_casadi_symbolic(casadi_symbolic_mode,['s_sot_' num2str(k)],1);
            w = {w{:}, s_sot_k};
            % intialize speed of time (sot), lower and upper bounds
            w0 = [w0; s_sot0];
            ubw = [ubw;s_sot_max];
            lbw = [lbw;s_sot_min];
            % index colector for sot variables
            ind_sot = [ind_sot,ind_total(end)+1:ind_total(end)+1];
            ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];
            J_regularize_sot = J_regularize_sot+(s_sot_k-1)^2;
        else
            if k == 0
                % only once
                s_sot_k = define_casadi_symbolic(casadi_symbolic_mode,['s_sot_' num2str(k+1)],1);
                w = {w{:}, s_sot_k};
                % intialize speed of time (sot), lower and upper bounds
                w0 = [w0; s_sot0];
                lbw = [lbw;s_sot_min];
                ubw = [ubw;s_sot_max];
                % index colector for sot variables
                ind_sot = [ind_sot,ind_total(end)+1:ind_total(end)+1];
                ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];
                % enforce gradient steps towward smaller values of s_sot to aid convergence (large variaties of s_sot = high nonlinearity)
                J_regularize_sot = J_regularize_sot+(s_sot_k)^2;
            end
        end
    else
        s_sot_k = 1; % speed of time is one, if no transfomrations used
    end

    %% General Nonlinear constraint (on control interval boundary)
    % The CasADi function g_ineq_fun and its lower and upper bound are provieded in model.
    if g_ineq_constraint
        g_ineq_k = g_ineq_fun(X_ki,Uk);
        g = {g{:}, g_ineq_k};
        lbg = [lbg; g_ineq_lb];
        ubg = [ubg; g_ineq_ub];
    end
    % path complementarity constraints
    if g_comp_path_constraint
        g_comp_path_k = g_comp_path_fun(X_ki,Uk)-sigma_p;
        g = {g{:}, g_comp_path_k};
        lbg = [lbg; g_comp_path_lb];
        ubg = [ubg; g_comp_path_ub];
    end

    sum_h_ki_control_interval_k = 0; % Integral of the clock state (if not time freezing) on the current stage.

    %% Least square terms
    J = J+T/N_stages*f_lsq_x_fun(X_ki,x_ref_val(:,k+1));
    if n_u > 0
        J = J+T/N_stages*f_lsq_u_fun(Uk,u_ref_val(:,k+1));
    end
    %% Loop over all finite elements in the current k-th control stage.
    for i = 0:N_finite_elements(k+1)-1
        %%  Sum of lambda and theta for current finite elememnt
        sum_Theta_ki = 0;  % initialize sum of theta's (the pint at t_n is not included)
        sum_Lambda_ki = Lambda_end_previous_fe;
        %% Step size in FESD, Speed of Time variables, Step equilibration constraints
        if ~use_fesd
            h_ki = h_k(k+1);
            if time_optimal_problem && ~use_speed_of_time_variables
                % final time is free, which implies a free step-size
                h_ki = T_final/(N_stages*N_finite_elements(k+1));
            end
        else
            % Define step-size variables, if FESD is used.
            if  k>0 || i>0
                h_ki_previous = h_ki;
            end
            % define step-size
            h_ki = define_casadi_symbolic(casadi_symbolic_mode,['h_'  num2str(k) '_' num2str(i)],1);
            w = {w{:}, h_ki};
            w0 = [w0; h0_k(k+1)];
            ubw = [ubw;ubh(k+1)];
            lbw = [lbw;lbh(k+1)];
            % index sets for step-size variables
            ind_h = [ind_h,ind_total(end)+1:ind_total(end)+1];
            ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];


            %             if (k > 0 && (i > 0 || couple_across_stages))
            if i > 0
                delta_h_ki = h_ki - h_ki_previous;
            else
                delta_h_ki  = 0;
            end

            % Terms for heuristic step equilibration
            if heuristic_step_equilibration
                switch heuristic_step_equilibration_mode
                    case 1
                        J_regularize_h  = J_regularize_h + (h_ki-h_k(k+1))^2;
                    case 2
                        J_regularize_h  = J_regularize_h + delta_h_ki^2;
                    otherwise
                        error('Pick heuristic_step_equilibration_mode between 1 and 2.');
                end
            end
        end

        integral_clock_state_k = integral_clock_state_k + h_ki*s_sot_k;
        sum_h_ki_control_interval_k = sum_h_ki_control_interval_k + h_ki;
        sum_h_ki_all = sum_h_ki_all + h_ki;

        %% Define Variables at stage points (IRK stages) for the current finite elements
        switch irk_representation
            case 'integral'
                X_ki_stages = {};
            case 'differential'
                V_ki_stages = {}; % for variables for derivative stage values
                X_ki_stages = {}; % for symbolic expressions of state stage values
        end
        %algebraic states
        Z_ki_stages = {}; % collects the vector of x and z at every irk stage

        % Collect in this struct the current algebraic variables, they are needed for cross complementarities
        Lambda_ki = {};
        Theta_ki = {};
        Mu_ki = {};

        % loop over all stage points to carry out defintions, initializations and bounds
        for j=1:n_s
            switch irk_representation
                case 'integral'
                    % define symbolic variables for values of diff. state a stage points
                    X_ki_stages{j} = define_casadi_symbolic(casadi_symbolic_mode,['X_'  num2str(k) '_' num2str(i) '_' num2str(j) ],n_x);
                    w = {w{:}, X_ki_stages{j}};
                    w0 = [w0; x0];
                    ind_x = [ind_x,ind_total(end)+1:ind_total(end)+n_x];
                    ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_x];
                    if x_box_at_stg
                        lbw = [lbw; lbx];
                        ubw = [ubw; ubx];
                    else
                        lbw = [lbw; -inf*ones(n_x,1)];
                        ubw = [ubw; inf*ones(n_x,1)];
                    end
                case 'differential'
                    V_ki_stages{j} = define_casadi_symbolic(casadi_symbolic_mode,['V_'  num2str(k) '_' num2str(i) '_' num2str(j) ],n_x);
                    w = {w{:}, V_ki_stages{j}};
                    w0 = [w0; v0];
                    ind_v = [ind_v,ind_total(end)+1:ind_total(end)+n_x];
                    ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_x];
                    lbw = [lbw; -inf*ones(n_x,1)];
                    ubw = [ubw; inf*ones(n_x,1)];

                    if lift_irk_differential
                        X_ki_stages{j} = define_casadi_symbolic(casadi_symbolic_mode,['X_'  num2str(k) '_' num2str(i) '_' num2str(j)],n_x);

                        w = {w{:}, X_ki_stages{j}};
                        w0 = [w0; x0];
                        ind_x = [ind_x,ind_total(end)+1:ind_total(end)+n_x];
                        ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_x];
                        % Bounds
                        if x_box_at_stg
                            lbw = [lbw; lbx];
                            ubw = [ubw; ubx];
                        else
                            lbw = [lbw; -inf*ones(n_x,1)];
                            ubw = [ubw; +inf*ones(n_x,1)];
                        end
                    end
            end
            % Note that the algebraic variablies are treated the same way in both irk representation modes.
            if strcmp(casadi_symbolic_mode, 'SX') || strcmp(casadi_symbolic_mode, 'casadi.SX') % TODO: remove this or
                zkij = [];
                for iz =1:n_z
                    zkij = vertcat(zkij, define_casadi_symbolic(casadi_symbolic_mode,[name(model.z(iz)) '_' num2str(j)], 1));
                end
                Z_ki_stages{j} = zkij;
            else
                Z_ki_stages{j} = define_casadi_symbolic(casadi_symbolic_mode,['Z_'  num2str(k) '_' num2str(i) '_' num2str(j)],n_z);
            end

            w = {w{:}, Z_ki_stages{j}};
            % index sets
            ind_z = [ind_z,ind_total(end)+1:ind_total(end)+n_z];
            ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_z];
            % bounds and initialization
            lbw = [lbw; lbz];
            ubw = [ubw; ubz];
            w0 = [w0; z0];

            % collection of all lambda and theta for current finite element, they are used for cross complementarites and step equilibration
            switch pss_mode
                case 'Stewart'
                    Theta_kij = Z_ki_stages{j}(1:n_theta);
                    Lambda_kij = Z_ki_stages{j}(n_theta+1:2*n_theta);
                    Mu_kij = Z_ki_stages{j}(2*n_theta+1);
                case 'Step'
                    Theta_kij = [Z_ki_stages{j}(1:n_alpha);e_alpha-Z_ki_stages{j}(1:n_alpha)];
                    Lambda_kij = Z_ki_stages{j}(n_alpha+1:3*n_alpha);
                    Mu_kij = [];
            end
            Theta_ki = {Theta_ki{:}, Theta_kij};
            Lambda_ki = {Lambda_ki{:}, Lambda_kij};
            Mu_ki = {Mu_ki{:}, Mu_kij};
            % Sum \theta and \lambda over the current finite element.
            if use_fesd
                sum_Lambda_ki =  sum_Lambda_ki + Lambda_kij;
                sum_Theta_ki =  sum_Theta_ki  +  Theta_kij;
            end
            % Update the standard complementarity
            J_comp_std = J_comp_std + J_cc_fun(Z_ki_stages{j});
        end
        % differential stage values (either as sym expresions or lifted variables)
        % define symbolic expression for diff state values at stage points. note that they are not degrees of
        % freedom (V_ki_stages) are, they are just used to increase readiblitiy.
        % X_{k,n,j} = x_n + h_n \sum_{i=1}^{n_s} a_{j,i}v_{n,j}
        if isequal(irk_representation,'differential')
            for j = 1:n_s
                x_temp = X_ki;
                % irk equations for stages
                for r = 1:n_s
                    x_temp = x_temp + h_ki*A_irk(j,r)*V_ki_stages{r};
                end
                if lift_irk_differential
                    X_ki_lift{j} =  X_ki_stages{j} - x_temp;
                else
                    X_ki_stages{j} = x_temp;
                end
            end
        end

        %% Additional variables in case of schemes not containting the boundary point as stage point, e.g., Gauss-Legendre schemes
        if use_fesd && ~right_boundary_point_explicit &&  (k<N_stages-1 || i< N_finite_elements(k+1)-1)
            switch pss_mode
                case 'Stewart'
                    Lambda_ki_end = define_casadi_symbolic(casadi_symbolic_mode,['Lambda_' num2str(k) '_' num2str(i) '_end'],n_theta);
                    Mu_ki_end = define_casadi_symbolic(casadi_symbolic_mode,['Mu_' num2str(k) '_' num2str(i) '_end'],n_simplex);
                    w = {w{:}, Lambda_ki_end};
                    w = {w{:}, Mu_ki_end};
                    % bounds and index sets
                    w0 = [w0; z0(n_theta+1:end)];
                    lbw = [lbw; lbz(n_theta+1:end)];
                    ubw = [ubw; ubz(n_theta+1:end)];
                    ind_boundary = [ind_boundary,ind_total(end)+1:(ind_total(end)+n_z-n_theta)];
                    ind_total  = [ind_total,ind_total(end)+1:(ind_total(end)+n_z-n_theta)];
                    Lambda_ki = {Lambda_ki{:}, Lambda_ki_end};
                    Mu_ki = {Mu_ki{:}, Mu_ki_end};
                case 'Step'
                    Lambda_ki_end = define_casadi_symbolic(casadi_symbolic_mode,['Lambda_' num2str(k) '_' num2str(i) '_end'],2*n_alpha);
                    Mu_ki_end = [];
                    w = {w{:}, Lambda_ki_end};
                    % bounds and index sets
                    w0 = [w0; z0(n_alpha+1:3*n_alpha)];
                    lbw = [lbw; lbz(n_alpha+1:3*n_alpha)];
                    ubw = [ubw; ubz(n_alpha+1:3*n_alpha)];
                    ind_boundary = [ind_boundary,ind_total(end)+1:(ind_total(end)+2*n_alpha)];
                    ind_total  = [ind_total,ind_total(end)+1:(ind_total(end)+2*n_alpha)];
                    Lambda_ki = {Lambda_ki{:}, Lambda_ki_end};
                    Mu_ki = {Mu_ki{:}, []};
            end
            % For cross comp
            sum_Lambda_ki =  sum_Lambda_ki + Lambda_ki_end;
        end

        %% Update struct with all complementarity related quantities
        if use_fesd
            comp_var_current_fe.sum_Lambda_ki = sum_Lambda_ki;
            comp_var_current_fe.sum_Theta_ki= sum_Theta_ki;
            comp_var_current_fe.Lambda_end_previous_fe = Lambda_end_previous_fe;
        end
        comp_var_current_fe.Lambda_ki = Lambda_ki;
        comp_var_current_fe.Theta_ki = Theta_ki;

        %% Continiuity of lambda, the boundary values of lambda and mu
        if use_fesd
            if right_boundary_point_explicit
                Z_kd_end = Z_ki_stages{n_s};
                switch pss_mode
                    case 'Stewart'
                        Lambda_ki_end = Z_kd_end(n_theta+1:2*n_theta);
                    case 'Step'
                        Lambda_ki_end = Z_kd_end(n_alpha+1:3*n_alpha);
                end
            else
                switch pss_mode
                    case 'Stewart'
                        Z_kd_end = [zeros(n_theta,1);Lambda_ki_end;Mu_ki_end];
                    case 'Step'
                        Z_kd_end = [zeros(n_alpha,1);Lambda_ki_end;Mu_ki_end;zeros(n_beta,1);zeros(n_gamma,1)];
                end

            end
            % Update lambda previous at finite element level
            Lambda_end_previous_fe = Lambda_ki_end;
        end

        %% The IRK Equations: evaluate equations (dynamics, algebraic, complementarities standard and cross at every stage point)
        % initialization with value of at at left boundary point
        switch irk_representation
            case 'integral'
                Xk_end = D(1)*X_ki;
                % Note that the polynomial is initialized with the previous value % (continuity connection)
                % Xk_end = D(1)*Xk;   % X_k+1 = X_k0 + sum_{j=1}^{d} D(j)*X_kj; Xk0 = D(1)*Xk
                % X_k+1 = X_k0 + sum_{j=1}^{n_s} D(j)*X_kj; Xk0 = D(1)*Xk
            case 'differential'
                Xk_end = X_ki; % initialize with x_n;
        end

        for j=1:n_s
            switch irk_representation
                case 'integral'
                    % Expression for the state derivative at the stage point
                    xp = C(1,j+1)*X_ki;
                    % Lagrange polynomial with values at state
                    for r=1:n_s
                        xp = xp + C(r+1,j+1)*X_ki_stages{r};
                    end
                    % Evaluate Differential and Algebraic Equations at stage points
                    if n_u > 0
                        [fj, qj] = f_x_fun(X_ki_stages{j},Z_ki_stages{j},Uk);
                        gj = g_z_all_fun(X_ki_stages{j},Z_ki_stages{j},Uk);
                    else
                        [fj, qj] = f_x_fun(X_ki_stages{j},Z_ki_stages{j});
                        gj = g_z_all_fun(X_ki_stages{j},Z_ki_stages{j});
                    end
                    if time_rescaling && use_speed_of_time_variables
                        % rescale equations
                        fj = s_sot_k*fj;
                        qj = s_sot_k*qj;
                        gj = s_sot_k*gj;
                    end
                    % Add contribution to the end state, Attention qj changes with time scaling!
                    Xk_end = Xk_end + D(j+1)*X_ki_stages{j};
                    % Add contribution to quadrature function
                    J = J + B(j+1)*qj*h_ki;
                    % Append IRK equations to NLP constraint
                    g = {g{:}, h_ki*fj - xp};

                case 'differential'
                    % Evaluate Differential and Algebraic Equations at stage points
                    if n_u > 0
                        [fj, qj] = f_x_fun(X_ki_stages{j},Z_ki_stages{j},Uk);
                        gj = g_z_all_fun(X_ki_stages{j},Z_ki_stages{j},Uk);
                    else
                        [fj, qj] = f_x_fun(X_ki_stages{j},Z_ki_stages{j});
                        gj = g_z_all_fun(X_ki_stages{j},Z_ki_stages{j});
                    end

                    if time_rescaling && use_speed_of_time_variables
                        % rescale equations
                        fj = s_sot_k*fj;
                        qj = s_sot_k*qj;
                        gj = s_sot_k*gj;
                    end
                    % Add contribution to the end state and quadrature term
                    if use_fesd
                        Xk_end = Xk_end + h_ki*b_irk(j)*V_ki_stages{j};
                        J = J + h_ki*b_irk(j)*qj;
                    else
                        Xk_end = Xk_end + h_k(k+1)*b_irk(j)*V_ki_stages{j};
                        J = J + h_k(k+1)*b_irk(j)*qj;
                    end
                    % Append IRK equations (differential part) to NLP constraint, note that we dont distingiush if use_fesd is on/off
                    % since it was done in the defintion of X_ki_stages which enters f_j
                    g = {g{:}, fj - V_ki_stages{j}};

                    % lifting considerations
                    if lift_irk_differential
                        g = {g{:}, X_ki_lift{j}};
                        lbg = [lbg; zeros(n_x,1)];
                        ubg = [ubg; zeros(n_x,1)];
                    end
                    if  x_box_at_stg && ~lift_irk_differential
                        g = {g{:};X_ki_stages{j}};
                        lbg = [lbg; lbx];
                        ubg = [ubg; ubx];
                    end
            end
            % Append IRK equations (algebraic part) to NLP constraint (same for both representations)
            g = {g{:}, gj};
            lbg = [lbg; zeros(n_x,1); zeros(n_algebraic_constraints,1)];
            ubg = [ubg; zeros(n_x,1); zeros(n_algebraic_constraints,1)];

            %% General nonlinear constraint at stage points
            if g_ineq_constraint && g_ineq_at_stg
                % indepednet of the fact is it lifter od not in the
                % differential case
                g_ineq_k = g_ineq_fun(X_ki_stages{j},Uk);
                g = {g{:}, g_ineq_k};
                lbg = [lbg; g_ineq_lb];
                ubg = [ubg; g_ineq_ub];
            end

            % path complementarity constraints
            if g_comp_path_constraint   && g_ineq_at_stg
                g_comp_path_k = g_comp_path_fun(X_ki_stages{j},Uk)-p(1);
                g = {g{:}, g_comp_path_k};
                lbg = [lbg; g_comp_path_lb];
                ubg = [ubg; g_comp_path_ub];
            end

            %% Complementarity constraints (standard and cross)
            % Prepare Input for Cross Comp Function
            n_cross_comp_i = 0;
            % Update current index
            current_index.k = k;  current_index.i = i; current_index.j = j;
            % Create cross comp constraints
            results_cross_comp = create_complementarity_constraints(use_fesd,cross_comp_mode,comp_var_current_fe,dimensions,current_index);
            g_cross_comp_j  = results_cross_comp.g_cross_comp_j;
            % Store updatet sums

            if i == N_finite_elements(k+1)-1 && j == n_s
                % Restart sum at end of the current control interval
                comp_var_current_fe.cross_comp_k = 0;
            else
                comp_var_current_fe.cross_comp_k = results_cross_comp.cross_comp_k;
            end
            comp_var_current_fe.cross_comp_all = results_cross_comp.cross_comp_all;
            % Dimensions of cross/std complementarities
            n_cross_comp_j = length(g_cross_comp_j);
            n_cross_comp_i = n_cross_comp_i+n_cross_comp_j;
            n_cross_comp(i+1,k+1) = n_cross_comp_i;
            g_all_comp_j = [g_cross_comp_j];
            n_all_comp_j = length(g_all_comp_j);

            %% Reformulation/relaxation of complementarity constraints
            % Treatment and reformulation of all complementarity constraints (standard and cross complementarity), their treatment depends on the chosen MPCC Method.
            mpcc_var_current_fe.J = J;
            mpcc_var_current_fe.g_all_comp_j = g_all_comp_j;
            if s_ell_inf_elastic_exists
                mpcc_var_current_fe.s_elastic = s_elastic;
            end
            if s_ell_1_elastic_exists
                s_elastic_ell_1  = define_casadi_symbolic(casadi_symbolic_mode,['s_elastic_'  num2str(k) '_' num2str(i) '_' num2str(j)], n_all_comp_j);
                w = {w{:}, s_elastic_ell_1};
                ind_elastic = [ind_elastic,ind_total(end)+1:ind_total(end)+n_all_comp_j];
                ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_all_comp_j];
                lbw = [lbw; s_elastic_min*ones(n_all_comp_j,1)];
                ubw = [ubw; s_elastic_max*ones(n_all_comp_j,1)];
                w0 = [w0;s_elastic_0*ones(n_all_comp_j,1)];
                % sum of all elastic variables, to be passed penalized
                sum_s_elastic = sum_s_elastic+ sum(s_elastic_ell_1);
                % pass to constraint creation
                mpcc_var_current_fe.s_elastic = s_elastic_ell_1;
            end

            [J,g_comp,g_comp_lb,g_comp_ub] = reformulate_mpcc_constraints(objective_scaling_direct,mpcc_mode,mpcc_var_current_fe,dimensions,current_index);
            g = {g{:}, g_comp};
            lbg = [lbg; g_comp_lb];
            ubg = [ubg; g_comp_ub];
        end

        %% Step equilibration
        step_equilibration_constrains;

        %% Continuity condition - new NLP variable for state at end of a finite element
        % Conntinuity conditions for differential state
        X_ki = define_casadi_symbolic(casadi_symbolic_mode,['X_'  num2str(k+1) '_' num2str(i+1) ],n_x);
        w = {w{:}, X_ki};
        w0 = [w0; x0];

        ind_x = [ind_x,ind_total(end)+1:ind_total(end)+n_x];
        ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_x];

        % box constraint always done at control interval boundary
        if x_box_at_fe
            lbw = [lbw; lbx];
            ubw = [ubw; ubx];
        else
            lbw = [lbw; -inf*ones(n_x,1)];
            ubw = [ubw; -inf*ones(n_x,1)];
        end

        % Add equality constraint
        g = {g{:}, Xk_end-X_ki};
        lbg = [lbg; zeros(n_x,1)];
        ubg = [ubg; zeros(n_x,1)];

        %% Evaluate inequality constraints at finite elements boundaries
        % TODO?: This should be removed? the left boundary point is treated
        % after control defintion, the very last right point should be treated in the terminal constraint
        if g_ineq_constraint && g_ineq_at_fe && i<N_finite_elements(k+1)-1
            % the third flag is because at i = 0 the evaulation is at the control interval boundary (done above)
            g_ineq_k = g_ineq_fun(X_ki,Uk);
            g = {g{:}, g_ineq_k};
            lbg = [lbg; g_ineq_lb];
            ubg = [ubg; g_ineq_ub];
        end

        % path complementarity constraints
        if g_comp_path_constraint   && g_ineq_at_fe && i<N_finite_elements(k+1)-1
            g_comp_path_k = g_comp_path_fun(X_ki,Uk)-p(1);
            g = {g{:}, g_comp_path_k};
            lbg = [lbg; g_comp_path_lb];
            ubg = [ubg; g_comp_path_ub];
        end

        %% g_z_all constraint for boundary point and continuity of algebraic variables.
        if ~right_boundary_point_explicit && use_fesd && (k< N_stages-1 || i< N_finite_elements(k+1)-1)
            if n_u > 0
                temp = g_z_all_fun(X_ki,Z_kd_end,Uk);
            else
                temp = g_z_all_fun(X_ki,Z_kd_end);
            end
            gj = temp(1:end-n_lift_eq);
            lbg = [lbg; zeros(n_algebraic_constraints-n_lift_eq,1)];
            ubg = [ubg; zeros(n_algebraic_constraints-n_lift_eq,1)];
            g = {g{:}, gj};
        end
    end

    %% equidistant grid in numerical time (Multiple-shooting type discretization)
    if equidistant_control_grid && use_fesd
        if ~time_optimal_problem
            % numerical time is fixed.
            g = {g{:}, sum_h_ki_control_interval_k-h};
            lbg = [lbg; 0];
            ubg = [ubg; 0];
        else
            if ~time_freezing
                % if time_freezing on, time optimal problems are regarded via the clock state
                if use_speed_of_time_variables
                    % numerical time and speed of time are decoupled
                    g = {g{:}, sum_h_ki_control_interval_k-h};
                    lbg = [lbg; 0];
                    ubg = [ubg; 0];
                    g = {g{:}, integral_clock_state_k-T_final/N_stages};
                    lbg = [lbg; 0];
                    ubg = [ubg; 0];
                else
                    % numerical time and speed of time are lumped together
                    g = {g{:}, sum_h_ki_control_interval_k-T_final/N_stages};
                    lbg = [lbg; 0];
                    ubg = [ubg; 0];
                end
            end
        end
    end

    %% equidistant grid in phyisical time (Stage-wise constraints on the clock state)
    if time_freezing && stagewise_clock_constraint
        % This makes mostly sense if time freezin is on. Imposes the constraints t((k+1)*h) = (k+1)*h+t_0 , x0(end) = t_0.
        if time_optimal_problem
            g = {g{:}, Xk_end(end)-(k+1)*(T_final/N_stages)+x0(end)};
        else
            g = {g{:}, Xk_end(end)-(k+1)*h+x0(end)};
        end
        lbg = [lbg; 0];
        ubg = [ubg; 0];
        ind_g_clock_state = [ind_g_clock_state;length(lbg)];
    end
    %% Update integral of clock state in numerical time;
    integral_clock_state = integral_clock_state + integral_clock_state_k ; 
    integral_clock_state_k = 0;
end

%% Scalar-valued complementarity residual
if use_fesd
    % sum of all possible cross complementarities;
    J_comp_fesd = sum(results_cross_comp.cross_comp_all);
    J_comp =  J_comp_fesd;
else
    % no additional complementarites than the standard ones
    J_comp_fesd = J_comp_std;
    J_comp =  J_comp_std;
end

%%  -- Constraint for the terminal numerical and physical time (if no equidistant grids are required) --
% If the control grid is not equidistant, the constraint on sum of h happen only at the end.
% The constraints are splited to those which are related to numerical and physical time, to make it easier to read.

%% Terminal numerical and physical time
if time_freezing
    % Terminal Phyisical Time (Posssble terminal constraint on the clock state if time freezing is active).
    if time_optimal_problem
        g = {g{:}, Xk_end(end)-T_final};
        lbg = [lbg; 0];
        ubg = [ubg; 0];
    else
        if impose_terminal_phyisical_time && ~stagewise_clock_constraint
            g = {g{:}, Xk_end(end)-(T_ctrl_p)};
            lbg = [lbg; 0];
            ubg = [ubg; 0];
        else
            % no terminal constraint on the numerical time, as it
            % is implicitly determined by
        end
    end
    if equidistant_control_grid && ~stagewise_clock_constraint
        if ~time_optimal_problem
            g = {g{:}, Xk_end(end)-T};
            lbg = [lbg; 0];
            ubg = [ubg; 0];
        end
    end
else
    if ~use_fesd
        if time_optimal_problem
            % if time_freezing is on, everything is done via the clock state.
            if use_speed_of_time_variables
                g = {g{:}, integral_clock_state-T_final};
                lbg = [lbg; 0];
                ubg = [ubg; 0];
            else
                % otherwise treated via variable h_ki, i.e.,  h_ki =  T_final/(N_stages*N_FE)
            end
        end
    else
        % if equidistant_control_grid = true all time constraint are added in
        % the main control loop for every control stage k and the code
        % below is skipped
        if  ~equidistant_control_grid
            % T_num = T_phy = T_final =  T.
            % all step sizes add up to prescribed time T.
            % if use_speed_of_time_variables = true, numerical time is decupled from the sot scaling (no mather if local or not):
            if ~time_optimal_problem
                g = {g{:}, sum_h_ki_all-T};
                lbg = [lbg; 0];
                ubg = [ubg; 0];
            else
                if ~use_speed_of_time_variables
                    g = {g{:}, sum_h_ki_all-T_final};
                    lbg = [lbg; 0];
                    ubg = [ubg; 0];
                else
                    % T_num = T_phy = T_final \neq T.
                    g = {g{:}, sum_h_ki_all-T};
                    lbg = [lbg; 0];
                    ubg = [ubg; 0];
                    g = {g{:}, integral_clock_state-T_final};
                    lbg = [lbg; 0];
                    ubg = [ubg; 0];
                end
            end
        end
    end
end

%%  Time-optimal problems: the objective.
% Time-optimal problems: define auxilairy variable for the final time.
ind_t_final = [];
% if time_optimal_problem && (~use_speed_of_time_variables || time_freezing)
if time_optimal_problem
    % Add to the vector of unknowns
    w = {w{:}, T_final};
    w0 = [w0; T_final_guess];
    lbw = [lbw;T_final_min];
    ubw = [ubw;T_final_max];
    ind_t_final = [ind_total(end)+1];
    ind_total  = [ind_total,ind_total(end)+1];
    J = J + T_final;
end

%% Terminal least square term
J = J + f_lsq_T_fun(X_ki,x_ref_end_val);

%% Terminal Constraints
% Add Terminal Constraint
if terminal_constraint
    if relax_terminal_constraint_homotopy
        rho_terminal_p = 1/sigma_p;
    end

    g_terminal = g_terminal_fun(Xk_end);
    n_terminal = length(g_terminal_lb);
    if ~isequal(g_terminal_lb,g_terminal_lb)
        relax_terminal_constraint = 0;
        % inequality constraint are not relaxed
        if print_level >2
            fprintf('Info: Only terminal-equality constraint relaxation is supported, you have an inequality constraint.\n')
        end
    end
    switch relax_terminal_constraint
        case 0
            % hard constraint
            g = {g{:}, g_terminal };
            lbg = [lbg; g_terminal_lb];
            if relax_terminal_constraint_from_above
                ubg = [ubg; g_terminal_ub*0+inf];
            else
                ubg = [ubg; g_terminal_ub];
            end
        case 1
            % l_1
            % define slack variables
            s_terminal_ell_1= define_casadi_symbolic(casadi_symbolic_mode,'s_terminal_ell_1',n_terminal);
            w = {w{:}, s_terminal_ell_1};
            lbw = [lbw;-inf*ones(n_terminal,1)];
            ubw = [ubw;inf*ones(n_terminal,1)];
            w0 = [w0;1e3*ones(n_terminal,1)];
            ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_terminal];
            % relaxed constraints
            g = {g{:}, g_terminal-g_terminal_lb-s_terminal_ell_1};
            g = {g{:}, -(g_terminal-g_terminal_lb)-s_terminal_ell_1};
            lbg = [lbg; -inf*ones(2*n_terminal,1)];
            ubg = [ubg; zeros(2*n_terminal,1)];
            % penalize slack
            J = J + rho_terminal_p*sum(s_terminal_ell_1);
        case 2
            % l_2
            J = J + rho_terminal_p*(g_terminal-g_terminal_lb)'*(g_terminal-g_terminal_lb);
        case 3
            % l_inf
            % define slack variable
            s_terminal_ell_inf= define_casadi_symbolic(casadi_symbolic_mode,'s_terminal_ell_inf',1);
            w = {w{:}, s_terminal_ell_inf};
            lbw = [lbw;-inf];
            ubw = [ubw;inf];
            w0 = [w0;1e3];
            ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];
            % relaxed constraint
            g = {g{:}, g_terminal-g_terminal_lb-s_terminal_ell_inf*ones(n_terminal,1)};
            g = {g{:}, -(g_terminal-g_terminal_lb)-s_terminal_ell_inf*ones(n_terminal,1)};
            lbg = [lbg; -inf*ones(2*n_terminal,1)];
            ubg = [ubg; zeros(2*n_terminal,1)];
            % penalize slack
            J = J + rho_terminal_p*(s_terminal_ell_inf);
        case 4
            if exist('s_elastic','var')
                % l_inf
                % define slack variable
                % relaxed constraint
                g = {g{:}, g_terminal-g_terminal_lb-s_elastic*ones(n_terminal,1)};
                g = {g{:}, -(g_terminal-g_terminal_lb)-s_elastic*ones(n_terminal,1)};
                lbg = [lbg; -inf*ones(2*n_terminal,1)];
                ubg = [ubg; zeros(2*n_terminal,1)];
            else
                error('This mode of terminal constraint relaxation is only available if a MPCC elastic mode is used.')
            end
    end
end

% path complementarity constraints
if g_comp_path_constraint
    g_comp_path_k = g_comp_path_fun(Xk_end,Uk)-p(1);
    g = {g{:}, g_comp_path_k};
    lbg = [lbg; g_comp_path_lb];
    ubg = [ubg; g_comp_path_ub];
end

%% quadratic regularization for speed of time variables;
% should be off in time optimal problems and on in time-freezing.
if time_freezing
    J = J+rho_sot_p*J_regularize_sot;
end

%% Add Terminal Cost to objective
try
    J = J + f_q_T_fun(Xk_end);
catch
    warning('Terminal cost not defined');
end
J_objective = J;

%% Elastic mode variable for \ell_infty reformulations
if mpcc_mode >= 5 && mpcc_mode < 8
    % add elastic variable to the vector of unknowns and add objective contribution
    w = {w{:}, s_elastic};
    ind_elastic = [ind_elastic,ind_total(end)+1:ind_total(end)+1];
    ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];
    lbw = [lbw; s_elastic_min];
    ubw = [ubw; s_elastic_max];
    w0 = [w0;s_elastic_0];

    if objective_scaling_direct
        J = J + (1/sigma_p)*s_elastic;
    else
        J = sigma_p*J + s_elastic;
    end
end

%% barrier controled penalty formulation
if mpcc_mode >= 8 && mpcc_mode <= 10
    rho_elastic = define_casadi_symbolic(casadi_symbolic_mode,'rho_elastic',1);

    w = {w{:}, rho_elastic};
    ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];

    s_elastic_min = 1e-16;
    s_elastic_max = inf;
    rho_0 = max(rho_min,0.5);

    if convex_sigma_rho_constraint
        sigma_0 = (1+sigma_scale)*rho_scale*exp(-rho_lambda*rho_0);
    else
        sigma_0 = sigma_scale*rho_scale*exp(-rho_lambda*rho_0);
    end
    rho_max = (log(rho_scale)-log(s_elastic_min))/rho_lambda;

    lbw = [lbw; rho_min];
    ubw = [ubw; rho_max];
    w0 = [w0;rho_0];
    if nonlinear_sigma_rho_constraint
        if convex_sigma_rho_constraint
            g = {g{:}, -rho_scale*exp(-rho_lambda*rho_elastic)+s_elastic};
        else
            g = {g{:}, rho_scale*exp(-rho_lambda*rho_elastic)-s_elastic};
        end
    else
        g = {g{:}, -(rho_elastic-rho_max)-s_elastic};
    end

    lbg = [lbg; 0];
    ubg = [ubg; inf];

    % add elastic variable to the vector of unknowns and add objective contribution
    w = {w{:}, s_elastic};
    ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];
    lbw = [lbw; s_elastic_min];
    ubw = [ubw; s_elastic_max];
    w0 = [w0;s_elastic_0];
    J = J-rho_penalty*(rho_elastic^2)+sigma_penalty*s_elastic;
end
%% Elastic mode variable for \ell_1 reformulations
if mpcc_mode >= 11 && mpcc_mode <= 13
    if objective_scaling_direct
        J = J + (1/sigma_p)*sum_s_elastic;
    else
        J = sigma_p*J + sum_s_elastic;
    end
end
%% Objective Terms for Grid Regularization
% Heuristic Regularization.
if heuristic_step_equilibration || step_equilibration
    %     J = J + step_equilibration_penalty*J_regularize_h;
    J = J + rho_h_p*J_regularize_h;
end

%% CasADi Functions for objective complementarity residual
J_fun = Function('J_fun', {vertcat(w{:})},{J_objective});
comp_res = Function('comp_res',{vertcat(w{:}), p},{J_comp});
comp_res_fesd = Function('comp_res_fesd',{vertcat(w{:})},{J_comp_fesd});
comp_res_std = Function('comp_res_std',{vertcat(w{:})},{J_comp_std});

%% NLP Solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}),'p',p);
solver = nlpsol(solver_name, 'ipopt', prob,opts_ipopt);

%% Define CasADi function for the switch indicator function.
if step_equilibration
    nu_fun = Function('nu_fun', {vertcat(w{:})},{nu_vector});
end

%% settings update
settings.right_boundary_point_explicit  = right_boundary_point_explicit  ;
%% Outputs
model.prob = prob;
model.solver = solver;
model.g =  vertcat(g{:});
model.w =  vertcat(w{:});
model.p =  p;
model.J = J;
model.J_fun = J_fun;
model.comp_res = comp_res;
model.comp_res_fesd = comp_res_fesd;
model.comp_res_std = comp_res_std;

% create CasADi function for objective gradient.
nabla_J = J.jacobian(model.w);
nabla_J_fun = Function('nabla_J_fun', {vertcat(w{:}),p},{nabla_J});
model.nabla_J = nabla_J;
model.nabla_J_fun = nabla_J_fun;

if step_equilibration
    model.nu_fun = nu_fun;
end

% TODO: make member function
if print_level > 4
    disp("g")
    print_casadi_vector(vertcat(g{:}))
    disp('lbg, ubg')
    disp([lbg, ubg])

    disp("w")
    print_casadi_vector(vertcat(w{:}))
    disp('lbw, ubw')
    disp([lbw, ubw])
    disp('w0')
    disp(w0)

    disp('objective')
    disp(J)
end

%% Model update: all index sets and dimensions
model.ind_x = ind_x;
model.ind_elastic = ind_elastic;
model.ind_v = ind_v;
model.ind_z = ind_z;
model.ind_u = ind_u;
model.ind_h = ind_h;
model.ind_vf = ind_vf;
model.ind_g_clock_state = ind_g_clock_state;
model.ind_sot = ind_sot;
model.ind_t_final  = ind_t_final;
model.ind_boundary = ind_boundary;
model.n_cross_comp = n_cross_comp;
model.ind_total = ind_total;
model.h = h;
model.h_k = h_k;
model.p_val = p_val;
model.n_cross_comp_total = sum(n_cross_comp(:));
%% Store solver initialization data
solver_initialization.w0 = w0;
solver_initialization.lbw = lbw;
solver_initialization.ubw = ubw;
solver_initialization.lbg = lbg;
solver_initialization.ubg = ubg;
%% Output
varargout{1} = solver;
varargout{2} = solver_initialization;
varargout{3} = model;
varargout{4} = settings;

end
