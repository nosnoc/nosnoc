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
%% Experimental stuff
new_mpcc_treatment = 1;

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

%% Create mpcc functions for the bilinear term
% mpcc_function = 'bilinear';
mpcc_function = 'Fischer-Burmeister';
% mpcc_function = 'Chen-Chen-Kanzow';
mpcc_function = 'Natural-Residual';
% mpcc_function = 'Steffenson-Ulbrich';
% mpcc_function = 'Steffenson-Ulbrich-Pol';
mpcc_function ='Kanzow-Schwartz';
% mpcc_function = 'Lin-Fukushima';
% mpcc_function = 'Kadrani';

[Psi_mpcc_fun] = create_mpcc_function(mpcc_function,casadi_symbolic_mode);
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
    T_final0 = T;
end

%% Elastic Mode Variables
if mpcc_mode >= 5 && mpcc_mode <= 10
    S_elastic = define_casadi_symbolic(casadi_symbolic_mode,'s_elastic',1);
else
    S_elastic  = 0;
end

%% Formulate NLP - Start with an empty NLP
% degrees of freedom
w = [];
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
J_virtual_froces = 0;
% constraints
g = [];
lbg = [];
ubg = [];

X_k0 = define_casadi_symbolic(casadi_symbolic_mode,'X_0',n_x);
w = [w; X_k0];

if there_exist_free_x0
    x0_ub = x0;
    x0_lb = x0;
    x0_ub(ind_free_x0) = inf;
    x0_lb(ind_free_x0) = -inf;
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
ind_vf  = [];
ind_sot = []; % index for speed of time variable
ind_boundary = []; % index of bundary value lambda and mu
ind_total = ind_x;
% constraint index sets
ind_g_total = 0;
ind_g_mpcc = [];
ind_g_irk_x = [];
ind_g_irk_z = [];
ind_g_clock_state = [];

% Integral of the clock state if no time rescaling is present.
sum_h_ki = 0;
sum_h_ki_vec = [];
sum_h_ki_all = 0;

% vectors conting sum_lambda_ki sum_theta_ki
sum_Theta_ki_vec = [];
sum_Lambda_ki_vec = [];

% Initialization of forward and backward sums for the step-equilibration
sigma_theta_B_k = 0; % backward sum of theta at finite element k
sigma_lambda_B_k = 0; % backward sum of lambda at finite element k
sigma_lambda_F_k = 0; % forward sum of lambda at finite element k
sigma_theta_F_k = 0; % forward sum of theta at finite element k
nu_vector = [];
% Integral of the clock state if no time-freezing is present.
sum_S_sot =0;
n_cross_comp = zeros(max(N_finite_elements),N_stages);

% Continuity of lambda initialization
Z_kd_end = zeros(n_z,1);
Lambda_end_previous_fe = zeros(n_theta,1);

% initialize cross comp and mpcc related structs
mpcc_var_k.p = p;
comp_var_k.cross_comp_k = 0;
comp_var_k.cross_comp_all = 0;


% TODO: Polish this part
comp_var_k.g_cross_comp_k = zeros(n_theta,1);
comp_var_k.g_cross_comp_all = zeros(n_theta,1);
comp_var_k.cross_comp_residual_k = 0;
comp_var_k.cross_comp_residual_all = 0;

g_cross_comp_k = zeros(n_theta,1);
g_cross_comp_all = zeros(n_theta,1);
cross_comp_residual_k = 0;
cross_comp_residual_all = 0;



%% Formulate the NLP / Main Discretization loop
for k=0:N_stages-1
    %% Variables for the controls
    if n_u > 0
        U_k = define_casadi_symbolic(casadi_symbolic_mode,['U_' num2str(k)],n_u);
        w = [w; U_k];
        % intialize contros, lower and upper bounds
        w0 = [w0; u0];
        lbw = [lbw;lbu];
        ubw = [ubw;ubu];
        % index colector for contorl variables
        ind_u = [ind_u,ind_total(end)+1:ind_total(end)+n_u];
        ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_u];

        % add pseudo froces in time-freezing problesm to aid convergence
        if virtual_forces
            if penalize_virtual_forces
                J_virtual_froces = J_virtual_froces+f_q_virtual_fun(U_k);
            end

            if virtual_forces_convex_combination
                J_virtual_froces = J_virtual_froces+psi_vf^2;
            end
            if tighthen_virtual_froces_bounds
                J_virtual_froces= J_virtual_froces+f_q_virtual_fun(U_k);
                U_virtual_k = U_k(end-n_q+1:end);
                if mpcc_mode < 4
                    g = [g; U_virtual_k - M_virtual_forces*sigma_p];
                    g = [g; -M_virtual_forces*sigma_p-U_virtual_k ];
                else
                    g = [g;  U_virtual_k - M_virtual_forces*S_elastic];
                    g = [g;  -M_virtual_forces*S_elastic-U_virtual_k];
                end
                lbg = [lbg; -inf*ones(2*n_q,1)];
                ubg = [ubg; 0*ones(2*n_q,1)];
                ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+2*n_q];
            end
        end
    end

    %%  Time rescaling of the stages (speed of time) to acchieve e.g., a desired final time in Time-Freezing or to solve time optimal control problems.
    %     If only time_rescaling is true, then the h_k also make sure to addapt the length of the subintervals, if both
    %     time_rescaling && use_speed_of_time_variables are true, new control variables are introduecd, they can be per stage or one for the whole interval.
    if time_rescaling && use_speed_of_time_variables
        if k == 0 || local_speed_of_time_variable
            % at every stage
            S_sot_k = define_casadi_symbolic(casadi_symbolic_mode,['S_sot_' num2str(k)],1);
            w = [w; S_sot_k];
            % intialize speed of time (sot), lower and upper bounds
            w0 = [w0; s_sot0];
            ubw = [ubw;s_sot_max];
            lbw = [lbw;s_sot_min];
            % index colector for sot variables
            ind_sot = [ind_sot,ind_total(end)+1:ind_total(end)+1];
            ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];
            J_regularize_sot = J_regularize_sot+(S_sot_k-1)^2;
        end
    end

    %% General Nonlinear constraint (on control interval boundary)
    % The CasADi function g_ineq_fun and its lower and upper bound are provieded in model.
    if g_ineq_constraint
        g_ineq_k = g_ineq_fun(X_k0,U_k);
        g = [g;  g_ineq_k];
        lbg = [lbg; g_ineq_lb];
        ubg = [ubg; g_ineq_ub];
        ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+length(g_ineq_ub)];
    end

    %% Loop over all finite elements in the current k-th control stage.
    sum_h_ki = 0; % Integral of the clock state (if not time freezing) on the current stage.

    for i = 0:N_finite_elements(k+1)-1
        %%  Sum of lambda and theta for current finite elememnt
        sum_Theta_ki = 0;  % initalize sum of theta's (the pint at t_n is not included)
        sum_Lambda_ki = Lambda_end_previous_fe;
        %% Step size in FESD, Speed of Time variables, Step equilibration constraints
        if use_fesd
            % Define step-size variables, if FESD is used.
            if  k>0 || i>0
                h_ki_previous = h_ki;
            end
            % define step-size
            h_ki = define_casadi_symbolic(casadi_symbolic_mode,['h_'  num2str(k) '_' num2str(i)],1);
            w = [w; h_ki];
            w0 = [w0; h0_k(k+1)];
            ubw = [ubw;ubh(k+1)];
            lbw = [lbw;lbh(k+1)];
            sum_h_ki = sum_h_ki + h_ki;
            % index sets for step-size variables
            ind_h = [ind_h,ind_total(end)+1:ind_total(end)+1];
            ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];

            % integral of clock state (if no time-freezing is used, in time-freezing we use directly the (nonsmooth) differential clock state.
            if time_optimal_problem && use_speed_of_time_variables
                sum_S_sot = sum_S_sot + h_ki*S_sot_k;
            end
            % \Delta h_ki
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
                        error('Pick heuristic_step_equlibration_mode between 1 and 2.');
                end
            end
        else
            h_ki = h_k(k+1);
        end

        %% Define Variables at stage points (IRK stages) for the current finite elements
        switch irk_representation
            case 'integral'
                X_ki = define_casadi_symbolic(casadi_symbolic_mode,['X_'  num2str(k) '_' num2str(i)],[n_x,n_s]);
                w = [w; X_ki(:)];
                w0 = [w0; repmat(x0,n_s,1)];
                ind_x = [ind_x,ind_total(end)+1:ind_total(end)+n_x*n_s];
                ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_x*n_s];
                if x_box_at_stg
                    lbw = [lbw; repmat(lbx,n_s,1)];
                    ubw = [ubw; repmat(ubx,n_s,1)];
                else
                    lbw = [lbw; -inf*ones(n_s*n_x,1)];
                    ubw = [ubw; inf*ones(n_s*n_x,1)];
                end
            case 'differential'
                % derivatives at stage points
                V_ki = define_casadi_symbolic(casadi_symbolic_mode,['V_'  num2str(k) '_' num2str(i)],[n_x,n_s]);
                w = [w; V_ki(:)];
                w0 = [w0; repmat(v0,n_s,1)];
                lbw = [lbw; -inf*ones(n_s*n_x,1)];
                ubw = [ubw; inf*ones(n_s*n_x,1)];
                ind_v = [ind_v,ind_total(end)+1:ind_total(end)+n_x*n_s];
                ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_x*n_s];
                if lift_irk_differential
                    X_ki = define_casadi_symbolic(casadi_symbolic_mode,['X_'  num2str(k) '_' num2str(i)],[n_x,n_s]);
                    w = [w; X_ki(:)];
                    w0 = [w0; repmat(x0,n_s,1)];
                    ind_x = [ind_x,ind_total(end)+1:ind_total(end)+n_x*n_s];
                    ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_x*n_s];

                    if x_box_at_stg
                        lbw = [lbw; repmat(lbx,n_s,1)];
                        ubw = [ubw; repmat(ubx,n_s,1)];
                    else
                        lbw = [lbw; -inf*ones(n_s*n_x,1)];
                        ubw = [ubw; inf*ones(n_s*n_x,1)];
                    end
                    X_ki_lift_expressions = X_ki - (X_k0 + h_ki*V_ki*A_irk');
                else
                    % store in this matrix all expressions for state at stage points
                    X_ki = X_k0+h_ki*V_ki*A_irk';
                    if x_box_at_stg
                        g = [g; X_ki(:)];
                        lbg = [lbg; repmat(lbx,n_s,1) ];
                        ubg = [ubg; repmat(ubx,n_s,1) ];
                        ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+n_x*n_s];
                    end
                end
        end
        % algebraic variables
        Z_ki = define_casadi_symbolic(casadi_symbolic_mode,['Z_'  num2str(k) '_' num2str(i)],[n_z,n_s]);
        w = [w; Z_ki(:)];
        % index sets
        ind_z = [ind_z,ind_total(end)+1:ind_total(end)+n_z*n_s];
        ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_z*n_s];
        % bounds and initialization
        lbw = [lbw; repmat(lbz,n_s,1)];
        ubw = [ubw; repmat(ubz,n_s,1)];
        w0 = [w0; repmat(z0,n_s,1)];

        % collection of all lambda and theta for current finite element, they are used for cross complementarites and step equilibration
        switch pss_mode
            case 'Stewart'
                Theta_ki = Z_ki(1:n_theta,:);
                Lambda_ki = Z_ki(n_theta+1:2*n_theta,:);
                Mu_ki = Z_ki(2*n_theta+1,:);
            case 'Step'
                Theta_ki= [Z_ki(1:n_alpha,:);repmat(e_alpha,1,n_s)-Z_ki(1:n_alpha,:)];
                Lambda_ki = Z_ki(n_alpha+1:3*n_alpha,:);
                Mu_ki = [];
        end
        % Sum \theta and \lambda over the current finite element.
        sum_Lambda_ki =  sum_Lambda_ki + sum(Lambda_ki,2);
        sum_Theta_ki =  sum_Theta_ki + sum(Theta_ki,2);
        % Update the standard complementarity
        J_comp_std = J_comp_std + sum(J_cc_fun(Z_ki));

        %% Additional boundary points if c_{n_s} \neq 1. (cf. FESD paper)
        if use_fesd && (k<N_stages-1 || i< N_finite_elements(k+1)-1)
            if right_boundary_point_explicit
                switch pss_mode
                    case 'Stewart'
                        Lambda_ki_end = Z_ki(n_theta+1:2*n_theta,n_s);
                    case 'Step'
                        Lambda_ki_end = Z_ki(n_alpha+1:3*n_alpha,n_s);
                end
            else
                switch pss_mode
                    case 'Stewart'
                        % n_s + 1 -th variable within current finite element
                        Lambda_ki_end = define_casadi_symbolic(casadi_symbolic_mode,['Lambda_' num2str(k) '_' num2str(i) '_end'],n_theta);
                        Mu_ki_end = define_casadi_symbolic(casadi_symbolic_mode,['Mu_' num2str(k) '_' num2str(i) '_end'],n_simplex);
                        w = [w; Lambda_ki_end];
                        w = [w; Mu_ki_end];
                        % bounds and index sets
                        w0 = [w0; z0(n_theta+1:end)];
                        lbw = [lbw; lbz(n_theta+1:end)];
                        ubw = [ubw; ubz(n_theta+1:end)];
                        ind_boundary = [ind_boundary,ind_total(end)+1:(ind_total(end)+n_z-n_theta)];
                        ind_total  = [ind_total,ind_total(end)+1:(ind_total(end)+n_z-n_theta)];

                    case 'Step'
                        Lambda_ki_end = define_casadi_symbolic(casadi_symbolic_mode,['Lambda_' num2str(k) '_' num2str(i) '_end'],2*n_alpha);
                        Mu_ki_end = [];
                        w = [w; Lambda_ki_end];
                        % bounds and index sets
                        w0 = [w0; z0(n_alpha+1:3*n_alpha)];
                        lbw = [lbw; lbz(n_alpha+1:3*n_alpha)];
                        ubw = [ubw; ubz(n_alpha+1:3*n_alpha)];
                        ind_boundary = [ind_boundary,ind_total(end)+1:(ind_total(end)+2*n_alpha)];
                        ind_total  = [ind_total,ind_total(end)+1:(ind_total(end)+2*n_alpha)];
                end
                Lambda_ki = [Lambda_ki, Lambda_ki_end];
                Mu_ki= [Mu_ki, Mu_ki_end];
                sum_Lambda_ki = sum_Lambda_ki + Lambda_ki_end; % add boundary point for cross comp (Attention: Is it added twice?)
            end
        end
        % Defintion of all stage variables for current finite element is done.

        %% Update struct with all complementarity related quantities
        if use_fesd
            comp_var_k.sum_Lambda_ki = sum_Lambda_ki;
            comp_var_k.sum_Theta_ki= sum_Theta_ki;
            comp_var_k.Lambda_end_previous_fe = Lambda_end_previous_fe;
        end
        if new_mpcc_treatment
            if use_fesd 
                comp_var_k.Lambda_ki = [Lambda_end_previous_fe Lambda_ki];
            else 
                comp_var_k.Lambda_ki = [Lambda_ki];
            end
        else
            comp_var_k.Lambda_ki = Lambda_ki;
        end
        comp_var_k.Theta_ki = Theta_ki;

        %% Continiuity of lambda, the boundary values of lambda and mu  %% TODO: resolve its use for proper cross comp -- move there the Z_kde end
        if use_fesd
            % Update lambda previous at finite element level
            %             if i > 0 || k > 0
            Lambda_end_previous_fe = Lambda_ki_end;
            %             end
        end

        %% The IRK Equations: evaluate equations (dynamics, algebraic, complementarities standard and cross at every stage point)
        % Evaulate Differential and Algebraic Equations at stage points
        if n_u > 0
            if virtual_forces && virtual_forces_convex_combination
                if virtual_forces_parametric_multipler
                    [f_x, f_q] = f_x_fun(X_ki,Z_ki,U_k,sigma_p);
                else
                    [f_x, f_q] = f_x_fun(X_ki,Z_ki,U_k,psi_vf);
                end
            else
                [f_x, f_q] = f_x_fun(X_ki,Z_ki,U_k);
            end
            g_z_all = g_z_all_fun(X_ki,Z_ki,U_k);
        else
            [f_x, f_q] = f_x_fun(X_ki,Z_ki);
            g_z_all = g_z_all_fun(X_ki,Z_ki);
        end

        if time_rescaling && use_speed_of_time_variables
            % rescale equations
            f_x = S_sot_k*f_x;
            f_q = S_sot_k*f_q;
            g_z_all = S_sot_k*g_z_all;
        end

        switch irk_representation
            case 'integral'
                % All stage values in one matrix, i.e. interpolating points of collocation polynomial
                X_all = [X_k0 X_ki];
                % Get slope of interpolating polynomial (normalized)
                Pi_dot = X_all*C(:,2:end);
                % Match with ODE right-hand-side
                g_irk = Pi_dot - h_ki*f_x;
                J = J + h_ki*f_q*B(2:end);
                % State at end of finite elemnt
                Xk_end = X_all*D;
            case 'differential'
                g_irk = f_x - V_ki;
                J = J + h_ki*f_q*b_irk';
                % State at end of finite elemnt
                Xk_end = X_k0 + h_ki*V_ki*b_irk';
                if lift_irk_differential
                    g = [g;  X_ki_lift_expressions(:)];
                    lbg = [lbg; zeros(n_x*n_s,1)];
                    ubg = [ubg; zeros(n_x*n_s,1)];
                    ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+n_x*n_s];
                end
        end
        g = [g;  g_irk(:); g_z_all(:)];
        lbg = [lbg; zeros(n_x*n_s,1);zeros(n_algebraic_constraints*n_s,1)];
        ubg = [ubg; zeros(n_x*n_s,1);zeros(n_algebraic_constraints*n_s,1)];

        % Update index sets for RK equations
        ind_g_irk_x = [ind_g_irk_x,ind_g_total(end)+1:ind_g_total(end)+n_x*n_s];
        ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+n_x*n_s];

        ind_g_irk_z = [ind_g_irk_z,ind_g_total(end)+1:ind_g_total(end)+n_algebraic_constraints*n_s];
        ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+n_algebraic_constraints*n_s];

        %% General nonlinear constraint at stage points
        if g_ineq_constraint && g_ineq_at_stg
            % indepednet of the fact is it lifter od not in the
            % differential case
            g_ineq_k = g_ineq_fun(X_ki,U_k);
            g = [g;  g_ineq_k(:)];
            lbg = [lbg; repmat(g_ineq_lb,n_s,1)];
            ubg = [ubg; repmat(g_ineq_ub,n_s,1)];
            ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+length(g_ineq_lb)*n_s];
        end

        %% Complementarity constraints (standard and cross)
        for j=1:n_s
            % Prepare Input for Cross Comp Function
            n_cross_comp_i = 0;
            % Update current index
            current_index.k = k;  current_index.i = i; current_index.j = j;
            % Create cross comp constraints
            if ~new_mpcc_treatment
                results_cross_comp = create_complementarity_constraints(use_fesd,cross_comp_mode,comp_var_k,dimensions,current_index);
                g_cross_comp_j  = results_cross_comp.g_cross_comp_j;
                % Store updatet sums
                if i == N_finite_elements(k+1)-1 && j == n_s
                    % Restart sum at end of the current control interval
                    comp_var_k.cross_comp_k = 0;
                else
                    comp_var_k.cross_comp_k = results_cross_comp.cross_comp_k;
                end
                comp_var_k.cross_comp_all = results_cross_comp.cross_comp_all;
                % Dimensions of cross/std complementarities
                n_cross_comp_j = length(g_cross_comp_j);
                n_cross_comp_i = n_cross_comp_i+n_cross_comp_j;
                n_cross_comp(i+1,k+1) = n_cross_comp_i;
                g_all_comp_j = [g_cross_comp_j];
                n_all_comp_j = length(g_all_comp_j);
            end

            %% Reformulation/relaxation of complementarity constraints
            % Treatment and reformulation of all complementarity constraints (standard and cross complementarity), their treatment depends on the chosen MPCC Method.
            if ~new_mpcc_treatment
                mpcc_var_k.J = J;
                mpcc_var_k.g_all_comp_j = g_all_comp_j;
                mpcc_var_k.s_elastic = S_elastic;
                [J,g_comp,g_comp_lb,g_comp_ub] = reformulate_mpcc_constraints(objective_scaling_direct,mpcc_mode,mpcc_var_k,dimensions,current_index);
                g = [g;  g_comp];
                lbg = [lbg; g_comp_lb];
                ubg = [ubg; g_comp_ub];
                ind_g_mpcc = [ind_g_mpcc,ind_g_total(end)+1:ind_g_total(end)+n_cross_comp_j];
                ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+n_cross_comp_j];
            else
                %                 g = [g;  g_cross_comp_j];
                %                 lbg = [lbg; -inf*ones(n_cross_comp_j,1)];
                %                 ubg = [ubg; zeros(n_cross_comp_j,1)];
                %                 ind_g_mpcc = [ind_g_mpcc,ind_g_total(end)+1:ind_g_total(end)+n_cross_comp_j];
                %             ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+n_cross_comp_j];
            end

        end
        % New function for conmplementarity constraints (no loops and sums)
        if new_mpcc_treatment
            [results_cross_comp] = create_complementarity_constraints_rv(use_fesd,cross_comp_mode,comp_var_k,dimensions,Psi_mpcc_fun,sigma_p,current_index);
            unfold_struct(results_cross_comp,'caller');
                 if i == N_finite_elements(k+1)-1 
                    comp_var_k.g_cross_comp_k = zeros(n_theta,1);
                    comp_var_k.cross_comp_residual_k = 0;
                else
                    comp_var_k.cross_comp_residual_k = cross_comp_residual_k ;
                    comp_var_k.g_cross_comp_k = g_cross_comp_k;
                end
                comp_var_k.g_cross_comp_all = g_cross_comp_all; 
                comp_var_k.cross_comp_residual_all = cross_comp_residual_all;

            g_cross_comp_ki = results_cross_comp.g_cross_comp_ki;
            n_cross_comp_ki = results_cross_comp.n_cross_comp_ki;

            g = [g;  g_cross_comp_ki];
            lbg = [lbg; -inf*ones(n_cross_comp_ki,1)];
            ubg = [ubg; zeros(n_cross_comp_ki,1)];
            %                 % update index sets
            ind_g_mpcc = [ind_g_mpcc,ind_g_total(end)+1:ind_g_total(end)+n_cross_comp_ki];
            ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+n_cross_comp_ki];
        end
        %% Step equlibration
        step_equilibration_constrains;
        %% Continuity condition - new NLP variable for state at end of a finite element
        % Conntinuity conditions for differential state
        X_k0 = define_casadi_symbolic(casadi_symbolic_mode,['X_'  num2str(k+1) '_' num2str(i+1) ],n_x);
        w = [w; X_k0];
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
        g = [g;  Xk_end-X_k0];
        lbg = [lbg; zeros(n_x,1)];
        ubg = [ubg; zeros(n_x,1)];
        ind_g_irk_x = [ind_g_irk_x,ind_g_total(end)+1:ind_g_total(end)+n_x];
        ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+n_x];

        %% Evaluate inequality constraints at finite elements boundaries
        % TODO?: This should be removed? the left boundary point is treated
        % after control defintion, the very last right point should be treated in the terminal constraint
        if g_ineq_constraint && g_ineq_at_fe && i<N_finite_elements(k+1)-1
            % the third flag is because at i = 0 the evaulation is at the control interval boundary (done above)
            g_ineq_k = g_ineq_fun(X_k0,U_k);
            g = [g;  g_ineq_k];
            lbg = [lbg; g_ineq_lb];
            ubg = [ubg; g_ineq_ub];
            ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+length(g_ineq_lb)];
        end

        %% g_z_all constraint for boundary point and continuity of algebraic variables.
        if ~right_boundary_point_explicit && use_fesd && (k< N_stages-1 || i< N_finite_elements(k+1)-1)
            switch pss_mode
                case 'Stewart'
                    Z_kd_end = [zeros(n_theta,1);Lambda_ki_end;Mu_ki_end];
                case 'Step'
                    Z_kd_end = [zeros(n_alpha,1);Lambda_ki_end;Mu_ki_end;zeros(n_beta,1);zeros(n_gamma,1)];
            end
            if n_u > 0
                temp = g_z_all_fun(X_k0,Z_kd_end,U_k);
            else
                temp = g_z_all_fun(X_k0,Z_kd_end);
            end
            g_z_all = temp(1:end-n_lift_eq);
            lbg = [lbg; zeros(n_algebraic_constraints-n_lift_eq,1)];
            ubg = [ubg; zeros(n_algebraic_constraints-n_lift_eq,1)];
            g = [g;  g_z_all];
            ind_g_irk_z = [ind_g_irk_z,ind_g_total(end)+1:ind_g_total(end)+n_algebraic_constraints-n_lift_eq];
            ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+n_algebraic_constraints-n_lift_eq];
        end
    end
    sum_h_ki_vec = [sum_h_ki_vec,sum_h_ki];
    sum_Theta_ki_vec = [sum_Theta_ki_vec,sum_Theta_ki];
    sum_Lambda_ki_vec = [sum_Lambda_ki_vec,sum_Lambda_ki];
    %% equidistant grid in numerical time (Multiple-shooting type discretization)
    if equidistant_control_grid
        if time_rescaling && time_optimal_problem
            % the querry here is because: No Time freezing: time_opt => time_rescaling (so this is always true if time_rescaling is on)
            % If time freezing is present, but not ime optimal problem, then the final numerical time should not be changed, hence the query.
            if use_speed_of_time_variables
                % numerical time and speed of time are decoupled
                g = [g;  sum_h_ki-h];
                lbg = [lbg; 0];
                ubg = [ubg; 0];
                ind_g_total = [ind_g_total,ind_g_total(end)+1];
            else
                % numerical time and speed of time are lumped together
                g = [g;  sum_h_ki-T_final/N_stages];
                lbg = [lbg; 0];
                ubg = [ubg; 0];
                ind_g_total = [ind_g_total,ind_g_total(end)+1];
            end
        else
            % numerical time is fixed.
            g = [g;  sum_h_ki-h];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            ind_g_total = [ind_g_total,ind_g_total(end)+1];
        end
    end

    %% equidistant grid in phyisical time (Stage-wise constraints on the colock state)
    if time_freezing && stagewise_clock_constraint
        % This makes mostly sense if time freezin is on. Imposes the constraints t((k+1)*h) = (k+1)*h+t_0 , x0(end) = t_0.
        if time_optimal_problem
            % this seems to be a critical constraint, T_final has to be
            g = [g; Xk_end(end)-(k+1)*(T_final/N_stages)+x0(end)];
        else
            g = [g; Xk_end(end)-(k+1)*h+x0(end)];
        end
        lbg = [lbg; 0];
        ubg = [ubg; 0];
        ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+1];
        ind_g_clock_state = [ind_g_clock_state;length(lbg)];
    end
end
sum_h_ki_all = sum(sum_h_ki_vec);
ind_g_total(1)  = []; % drop the zero that was needed for initialization
%% Scalar-valued commplementarity residual
if use_fesd
    % sum of all possible cross complementarities;
%     J_comp_fesd = sum(results_cross_comp.cross_comp_all);
    J_comp_fesd = cross_comp_residual_all;
    J_comp =  J_comp_fesd;
else
    % no additional complementarites than the standard ones
    J_comp_fesd = J_comp_std;
    J_comp =  J_comp_std;
end

%%  -- Constraint for the terminal numerical and physical time (if no equidistant grids are required) --
% If the control grid is not equidistant, the constraint on sum of h happen only at the end.
% The constraints are splited to those which are related to numerical and physical time, to make it easier to read.

%% Terminal numerical Time
if  ~equidistant_control_grid
    if time_rescaling
        % this means this is a time optimal control problem without time freezing
        if use_speed_of_time_variables
            % numerical time is decupled from the sot scaling (no mather if local or not):
            g = [g; sum_h_ki_all-T];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            ind_g_total = [ind_g_total,ind_g_total(end)+1];
            % T = T_num,  T_phy = s_sot*T.
            % note that objective contribution is added below (as it happens no
            % mather if we use local or nonlocal variables and equdistiant grid or not.)
        else
            if time_optimal_problem
                % the querry here is because: No Time freezing: time_opt => time_rescaling (so this is always true if time_rescaling is on)
                % If time freezing is present, but not ime optimal problem, then the final numerical time should not be changed, hence the query.
                g = [g; sum_h_ki_all-T_final];
                lbg = [lbg; 0];
                ubg = [ubg; 0];
                ind_g_total = [ind_g_total,ind_g_total(end)+1];
                % add time to objective.
                J = J+T_final;
            end
            % T_num = T_final = T_phy \neq T.
        end
    else
        % fixed numerical time (no time freezing or time optimal control, T_phy = T_num = T)
        g = [g; sum_h_ki_all-T];
        lbg = [lbg; 0];
        ubg = [ubg; 0];
        ind_g_total = [ind_g_total,ind_g_total(end)+1];
    end
end

%% Terminal Phyisical Time (Posssble terminal constraint on the clock state if time freezing is active).
if time_freezing
    if time_optimal_problem
        g = [g; Xk_end(end)-T_final];
        lbg = [lbg; 0];
        ubg = [ubg; 0];
        ind_g_total = [ind_g_total,ind_g_total(end)+1];
    else
        if impose_terminal_phyisical_time && ~stagewise_clock_constraint
            g = [g; Xk_end(end)-(T_ctrl_p)];
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            ind_g_total = [ind_g_total,ind_g_total(end)+1];
        else
            % no terminal constraint on the numerical time, as it
            % is implicityl determined by
        end
    end
end

if equidistant_control_grid && ~stagewise_clock_constraint && time_freezing
    if ~time_optimal_problem
        g = [g; Xk_end(end)-T];
        lbg = [lbg; 0];
        ubg = [ubg; 0];
        ind_g_total = [ind_g_total,ind_g_total(end)+1];
    end
end

%%  Time-optimal problems: the objective.
if time_optimal_problem && use_speed_of_time_variables
    if time_freezing
        J = J+T_final;
    else
        % This is both variants local and global sot vars. Hence, it is written here together and not in the loops or queries above.
        % ( In this case,the trivial constraint T_final = sum_s_sot) is ommited and provided implicitly)
        g = [g; sum_S_sot-T_final];
        lbg = [lbg; 0];
        ubg = [ubg; 0];
        ind_g_total = [ind_g_total,ind_g_total(end)+1];
        J = J+T_final;
        %         J = J+sum_s_sot;
        % sum_s_sot is the integral of the clock state if no time freezing is present.
    end
end
% Time-optimal problems: define auxilairy variable for the final time.
ind_t_final = [];
% if time_optimal_problem && (~use_speed_of_time_variables || time_freezing)
if time_optimal_problem
    % Add to the vector of unknowns
    w = [w; T_final];
    w0 = [w0; T_final0];
    lbw = [lbw;T_final_min];
    ubw = [ubw;T_final_max];
    ind_t_final = [ind_total(end)+1];
    ind_total  = [ind_total,ind_total(end)+1];
    J = J + T_final;
end

%% Terminal Constraints
% Add Terminal Constrint
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
            g = [g; g_terminal];
            lbg = [lbg; g_terminal_lb];
            if relax_terminal_constraint_from_above
                ubg = [ubg; g_terminal_ub*0+inf];
            else
                ubg = [ubg; g_terminal_ub];
            end
            ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+n_terminal];
        case 1
            % l_1
            % define slack variables
            s_terminal_ell_1= define_casadi_symbolic(casadi_symbolic_mode,'s_terminal_ell_1',n_terminal);
            w = [w; s_terminal_ell_1];
            lbw = [lbw;-inf*ones(n_terminal,1)];
            ubw = [ubw;inf*ones(n_terminal,1)];
            w0 = [w0;1e3*ones(n_terminal,1)];
            ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_terminal];
            % relaxed constraints
            g = [g; g_terminal-g_terminal_lb-s_terminal_ell_1];
            g = [g; -(g_terminal-g_terminal_lb)-s_terminal_ell_1];
            lbg = [lbg; -inf*ones(2*n_terminal,1)];
            ubg = [ubg; zeros(2*n_terminal,1)];
            ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+2*n_terminal];
            % penalize slack
            J = J + rho_terminal_p*sum(s_terminal_ell_1);
        case 2
            % l_2
            J = J + rho_terminal_p*(g_terminal-g_terminal_lb)'*(g_terminal-g_terminal_lb);
        case 3
            % l_inf
            % define slack variable
            s_terminal_ell_inf= define_casadi_symbolic(casadi_symbolic_mode,'s_terminal_ell_inf',1);
            w = [w; s_terminal_ell_inf];
            lbw = [lbw;-inf];
            ubw = [ubw;inf];
            w0 = [w0;1e3];
            ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];
            % relaxed constraint
            g = [g; g_terminal-g_terminal_lb-s_terminal_ell_inf*ones(n_terminal,1)];
            g = [g; -(g_terminal-g_terminal_lb)-s_terminal_ell_inf*ones(n_terminal,1)];
            lbg = [lbg; -inf*ones(2*n_terminal,1)];
            ubg = [ubg; zeros(2*n_terminal,1)];
            ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+2*n_terminal];
            % penalize slack
            J = J + rho_terminal_p*(s_terminal_ell_inf);
        case 4
            if exist('s_elastic','var')
                % l_inf
                % define slack variable
                % relaxed constraint
                g = [g; g_terminal-g_terminal_lb-S_elastic*ones(n_terminal,1);-(g_terminal-g_terminal_lb)-S_elastic*ones(n_terminal,1)];
                lbg = [lbg; -inf*ones(2*n_terminal,1)];
                ubg = [ubg; zeros(2*n_terminal,1)];
                ind_g_total = [ind_g_total,ind_g_total(end)+1:ind_g_total(end)+2*n_terminal];
            else
                error('This mode of terminal constraint relaxation is only avilable if a MPCC elastic mode is used.')
            end
    end
end

%% quadratic regularization for speed of time variables;
% should be off in time optimal problems and on in time-freezing.
if time_freezing
    J = J+rho_sot_p*J_regularize_sot;
end

%% Virtual forces
if virtual_forces
    if objective_scaling_direct
        J = J + (1/sigma_p)*J_virtual_froces;
    else
        J = J+J_virtual_froces;
    end
    if virtual_forces_convex_combination
        % treat the convex multipler of
        w = [w; psi_vf];
        ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];
        ind_vf  = ind_total(end);
        w0 = [w0;1];
        if virtual_forces_parametric_multipler
            g = [g; psi_vf-sigma_p];
            lbg =[lbg;0];
            ubg =[ubg;0];
            ind_g_total = [ind_g_total,ind_g_total(end)+1];

            lbw = [lbw; -inf];
            ubw = [ubw; inf];
        else
            lbw = [lbw; 0];
            ubw = [ubw; 1];
        end
    end
end

%% Add Terminal Cost to objective
try
    J = J + f_q_T_fun(Xk_end);
catch
    warning('Terminal cost not defined');
end
J_objective = J;

%% Elastic mode variable for \ell_infty reformulations
if new_mpcc_treatment
    % MPCC Reformulation experimental
    nlp_data.J = J; nlp_data.g = g; nlp_data.lbg = lbg; nlp_data.ubg = ubg;
    nlp_data.w = w; nlp_data.w0= w0; nlp_data.lbw = lbw; nlp_data.ubw = ubw;
    nlp_data.sigma_p = sigma_p;

    ind_data.ind_g_total = ind_g_total; ind_data.ind_g_mpcc = ind_g_mpcc;

    [nlp_data] = mpcc_solution_method(nlp_data,ind_data,mpcc_mode);
    unfold_struct(nlp_data,'caller');    
%     lbw(ind_z) = -inf;
else
    if mpcc_mode >= 5 && mpcc_mode < 8
        % add elastic variable to the vector of unknowns and add objeective contribution
        w = [w; S_elastic];
        ind_elastic = [ind_elastic,ind_total(end)+1:ind_total(end)+1];
        ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];
        lbw = [lbw; s_elastic_min];
        ubw = [ubw; s_elastic_max];
        w0 = [w0;s_elastic_0];

        if objective_scaling_direct
            J = J + (1/sigma_p)*S_elastic;
        else
            J = sigma_p*J + S_elastic;
        end
    end
end
%% Objective Terms for Grid Regularization
% Heuristic Regularization.
if heuristic_step_equilibration || step_equilibration
    J = J + rho_h_p*J_regularize_h;
end

%% CasADi Functions for objective complementarity residual
J_fun = Function('J_fun', {w},{J_objective});
J_virtual_froces_fun = Function('J_virtual_froces_fun', {w},{J_virtual_froces});
comp_res = Function('comp_res',{w},{J_comp});
comp_res_fesd = Function('comp_res_fesd',{w},{J_comp_fesd});
comp_res_std = Function('comp_res_std',{w},{J_comp_std});

%% NLP Solver
prob = struct('f', J, 'x', w, 'g', g ,'p',p);
solver = nlpsol(solver_name, 'ipopt', prob,opts_ipopt);

%% Define CasADi function for the switch indicator function.
if step_equilibration
    nu_fun = Function('nu_fun',{w},{nu_vector});
end

%% settings update
settings.right_boundary_point_explicit  = right_boundary_point_explicit  ;
%% Outputs
model.prob = prob;
model.solver = solver;
model.g = g;
model.w = w;
model.J = J;
model.J_fun = J_fun;
model.J_virtual_froces_fun = J_virtual_froces_fun;
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
