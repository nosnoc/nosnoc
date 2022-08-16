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
function [varargout] = create_nlp_nosnoc_opti(varargin)
% Info on this function.
%% read data
model = varargin{1};
settings = varargin{2};
%% CasADi
import casadi.*
%% Reformulation of the PSS into a DCS
[settings] = refine_user_settings(settings);
[model,settings] = model_reformulation_nosnoc(model,settings);

%% Fillin missing settings with default settings
[settings] = fill_in_missing_settings(settings,model);

%% Load user settings and model details
unfold_struct(settings,'caller')
unfold_struct(model,'caller');
%% Bounds on step-size (if FESD active)
if use_fesd
    ubh = (1+gamma_h)*h_k;
    lbh = (1-gamma_h)*h_k;
    if time_rescaling && ~use_speed_of_time_variables
        % if only time_rescaling is true, speed of time and step size all lumped together, e.g., \hat{h}_{k,i} = s_n * h_{k,i}, hence the bounds need to be extended.
        ubh = (1+gamma_h)*h_k*s_sot_max;
        lbh = (1-gamma_h)*h_k/s_sot_min;
    end
    % initigal guess for the step-size
    h0_k = h_k.*ones(N_stages,1);
end

%% Butcher Tableu (differential and integral representation)
switch irk_representation
    case 'integral'
        [B,C,D,tau_root] = generatre_butcher_tableu_integral(n_s,irk_scheme);
        if tau_root(end) == 1
            right_boundary_point_explicit  = 1;
        else
            right_boundary_point_explicit  = 0;
        end
    case 'differential'
        [A_irk,b_irk,c_irk,order_irk] = generatre_butcher_tableu(n_s,irk_scheme);
        if c_irk(end) <= 1+1e-9 && c_irk(end) >= 1-1e-9
            right_boundary_point_explicit  = 1;
        else
            right_boundary_point_explicit  = 0;
        end
    otherwise
        error('Choose irk_representation either: ''integral'' or ''differential''')
end

%% Elastic Mode variables
% This will be defind inside the loops and functions for the treatment

%% Time optimal control (time transformations)
if time_optimal_problem
    % the final time in time optimal control problems
    T_final = define_casadi_symbolic(casadi_symbolic_mode,'T_final',1);
    T_final_guess = T;
end

%% Formulation of the NLP
opti = Opti();

%% Objective terms
J = 0;
J_comp = 0;
J_comp_std = 0;
J_comp_cross = 0;
J_comp_step_eq = 0;
J_regularize_h = 0;
J_regularize_sot = 0;

%% Inital value
X_ki = opti.variable(n_x);
if there_exist_free_x0
    x0_ub = x0;
    x0_lb = x0;
    x0_ub(ind_free_x0) = inf;
    x0_lb(ind_free_x0) = -inf;
    x0(ind_free_x0) = 0;
    opti.subject_to(x0_lb<= X_ki <=x0_ub);
else
    opti.subject_to(X_ki==x0);

end
opti.set_initial(X_ki, x0);

%% Index sets and collcet all variables into vectors
% (TODO: are these index sets needed at all in opti?)
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
% Collect all states/controls
%% TODO: Consider dropping the s in Xs, Us etc.
X_boundary = {X_ki}; % differentail on FE boundary
X = {X_ki}; % differential on all points
V = {}; % stage derivative values in differential mode
Z = {};  % algebraic - PSS
Z_DAE = {}; % algebraic - DAE
U = {}; % Controls
H = {}; % Step-sizes
S_sot = {}; % Speed-of time
S_elastic = {}; % Elastic

%% Other initalizations of sums
sum_h_ki = []; % a vector of N_ctrl x 1, collecting all sums
sum_h_ki_all = 0; %  %TODO: make this at the end as sum of the entries of sum_h_ki; Note that this is the integral of the clock state if no time-freezing is used.
clock_state_integral_smooth = 0; % clock_state_integral_smooth
% Initalization of forward and backward sums for the step-equilibration
sigma_theta_B_k = 0; % backward sum of theta at finite element k
sigma_theta_F_k = 0; % forward sum of theta at finite element k
sigma_lambda_B_k = 0; % backward sum of lambda at finite element k
sigma_lambda_F_k = 0; % forward sum of lambda at finite element k
nu_vector = [];

sum_lambda_all = [];
sum_theta_all = [];

% Continuity of lambda initalization
Lambda_end_previous_fe = zeros(n_theta,1);
% TODO: check does n_theta dimension make sense for 'Step'??
% What is this Z_kd_end??
Z_kd_end = zeros(n_z,1);
% % initialize cross comp and mpcc related structs
% mpcc_var_current_fe.p = p;
% comp_var_current_fe.cross_comp_control_interval_k = 0;
% comp_var_current_fe.cross_comp_control_interval_all = 0;
n_cross_comp = zeros(max(N_finite_elements),N_stages);
%% Index nomenclature
%  k - control interval  {0,...,N_ctrl-1}
%  i - finite element    {0,...,N_fe-1}
%  j - stage             {1,...,n_s}
% TODO rename all Xk to X_k, all Uk to U_k, etc. to be consistent and increase readibility.
%% Main NLP loop over control intevrals/stages
for k=0:N_ctrl-1
    %% Define discrete-time control varibles for control interval k
    if n_u >0
        U_k = opti.variable(n_u);
        U{end+1} = U_k;
        opti.subject_to(lbu <= U_k <=ubu);
        opti.set_initial(U_k, u0);
        ind_u = [ind_u,ind_total(end)+1:ind_total(end)+n_u];
        ind_total = [ind_total,ind_total(end)+1:ind_total(end)+n_u];
        if virtual_forces
            % Possibly implement the costs and constraints for virutal forces in time-freezing
        end
    end

    %% Time-transformation variables (either one for the whole horizon or one per control interval)
    if time_rescaling && use_speed_of_time_variables
        if k == 0 || local_speed_of_time_variable
            S_sot_k = opti.variable(n_u);
            S_sot{end+1} = {S_sot_k}; % Speed-of time
            opti.subject_to(s_sot_min <= S_sot_k <=s_sot_max);
            opti.set_initial(S_sot_k, s_sot0); %TODO: rename s_sot0 to S_sot0 in the approipate places
            % index colector for sot variables
            ind_sot = [ind_sot,ind_total(end)+1:ind_total(end)+1];
            ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];
            J_regularize_sot = J_regularize_sot+(S_sot_k-S_sot_nominal)^2;
        end
    end

    %% General nonlinear constriant evaluated at left boundary of the control interval
    if g_ineq_constraint
        g_ineq_k = g_ineq_fun(X_ki,U_k);
        opti.subject_to(g_ineq_lb <= g_ineq_k <=g_ineq_ub);
    end

    %% Loop over finite elements for the current control interval
    sum_h_ki_temp = 0; % initalize sum for current control interval
    %% TODO: answer the questions: does vectorizing this loop make sense as vectorizing the loop over the stages below?
    for i = 0:N_FE(k+1)-1
        %%  Sum of all theta and lambda for current finite elememnt
        % Note that these are vector valued functions, sum_lambda_ki
        % contains the right boundary point of the previous interval because of the continuity of lambda (this is essential for switch detection, cf. FESD paper)
        sum_Theta_ki = zeros(n_theta,1);
        sum_Lambda_ki = Lambda_end_previous_fe;
        %% Step-size variables - h_ki;
        if use_fesd
            % Define step-size variables, if FESD is used.
            if  k>0 || i>0
                h_ki_previous = h_ki;
            end
            h_ki = opti.variable(1);
            H{end+1} = h_ki ;
            opti.subject_to(lbh(k+1) <= h_ki <=ubh(k+1));
            opti.set_initial(h_ki, h0_k);
            % index sets for step-size variables
            ind_h = [ind_h,ind_total(end)+1:ind_total(end)+1];
            ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];
            sum_h_ki_temp = sum_h_ki_temp + h_ki;

            % delta_h_ki = h_ki -h_{k-1,i}. Obeserve that the delta_h are  computed for every control interval separatly
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
                        % TODO: in refine_settings heuristic_step_equlibration_mode >2, set it to heuristic_step_equlibration_mode = 1
                        error('Pick heuristic_step_equlibration_mode between 1 and 2.');
                end
            end
            % Integral of clock state (if time flows without stopping, i.e., time_freezing = 0)
            if time_optimal_problem && use_speed_of_time_variables
                % TODO: consider writhing this forumla vectorized after all loops
                clock_state_integral_smooth = clock_state_integral_smooth + h_ki*s_sot_k; % integral of clock state (if no time-freezing is used, in time-freezing we use directly the (nonsmooth) differential clock state.
            end
        end
        %% Differntial variables at stage points
        % Defintion of differential variables
        switch irk_representation
            case 'integral'
                X_ki_stages = opti.variable(n_x, n_s); % matrix with all stage values for current FE i in control interval k
                X{end+1} = X_ki_stages;
                if x_box_at_stg
                    opti.subject_to(lbx <= X_ki_stages <= ubx);
                end
                opti.set_initial(X_ki_stages, repmat(x0,1,n_s));
                % Index sets
                ind_x = [ind_x,ind_total(end)+1:ind_total(end)+n_x*n_s];
                ind_total = [ind_total,ind_total(end)+1:ind_total(end)+n_x*n_s];
            case 'differential'
                % derivatives at stage points
                V_ki_stages = opti.variable(n_x, n_s); % for variables for derivative stage values
                V{end+1} = V_ki_stages;
                opti.subject_to(-inf*ones(n_x,1) <= V_ki_stages <= inf*ones(n_x,1));
                opti.set_initial(V_ki_stages, repmat(v0,1,n_s));
                ind_v = [ind_v,ind_total(end)+1:ind_total(end)+n_x];
                ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_x];

                X_ki_lift_expressions = []; % for symbolic expressions if X_ki_stages are lifted.

                if lift_irk_differential
                    X_ki_stages = opti.variable(n_x, n_s); % matrix with all stage values for current FE i in control interval k
                    X{end+1} = X_ki_stages;
                    if x_box_at_stg
                        opti.subject_to(lbx <= X_ki_stages <= ubx);
                    end
                    opti.set_initial(X_ki_stages, repmat(x0,1,n_s));
                    % Index sets
                    ind_x = [ind_x,ind_total(end)+1:ind_total(end)+n_x*n_s];
                    ind_total = [ind_total,ind_total(end)+1:ind_total(end)+n_x*n_s];
                else
                    % store in this matrix all expressions for state at stage points
                    X_ki_stages = [];
                end

                % differentail state values at stage points - either via  symbolic variables or new degress of freedom (lifting is on)
                % expressions for the state variable values at stage points follow:
                    % X_{k,i,j} = x_ki0 + h_n \sum_{i=1}^{n_s} a_{j,i}v_{i,j}
                for j = 1:n_s
                  % inital value (left boundary point of current FE)
                    x_temp = X_ki;
                    for r = 1:n_s
                        if use_fesd
                            x_temp = x_temp + h_ki*A_irk(j,r)*V_ki_stages(:,r);
                        else
                            x_temp = x_temp + h_k(k+1)*A_irk(j,r)*V_ki_stages(:,r);
                        end
                    end
                    % TODO: run few examples to see does this work properly
                        if lift_irk_differential
                            X_ki_lift_expressions = [X_ki_lift_expressions;X_ki_stages(:,j) - x_temp];
                            % TODO: This is a matrix function, be careful when adding as constraibnt
                        else
                            X_ki_stages = [X_ki_stages,x_temp];
                        end               
                end
        end

        %% Defintion of algebraic variables
        % Note that the algebraic variablies are treated the same way in both irk representation modes.
        Z_ki = opti.variable(n_z, n_s);
        Z{end+1} = Z_ki;
        opti.subject_to(lbz <= Z_ki <=ubz);
        opti.set_initial(Z_ki, repmat(z0,1,n_s));
        % Index sets
        ind_z = [ind_z,ind_total(end)+1:ind_total(end)+n_z*n_s];
        ind_total = [ind_total,ind_total(end)+1:ind_total(end)+n_z*n_s];
        % collection of all lambda and theta for current finite element, they are used for cross complementarites and step equilibration
        % They are matrices, e.g. n_theta * n_s, that collect all stage values
        switch pss_mode
            case 'Stewart'
                Theta_ki = Z_ki(1:n_theta,:);
                Lambda_ki = Z_ki(n_theta+1:2*n_theta,:);
                Mu_ki = Z_ki(2*n_theta+1,:);
            case 'Step'
                Theta_ki = [Z_ki(1:n_alpha,:);e_alpha-Z_ki(1:n_alpha,:)];
                Lambda_ki = Z_ki(n_alpha+1:3*n_alpha,:);
                Mu_kij = [];
        end
        % Remark for myslef: Theta_ki_current_fe replaced with Theta_ki now
        
        % Sum \theta and \lambda over the current finite element.
        if use_fesd
            % TODO: Check is it properly summed in casadi symbolics with sum(\cdot,2)
            sum_Theta_ki = sum(Theta_ki,2);
            sum_Lambda_ki = Lambda_end_previous_fe+sum(Lambda_ki,2);
            % these sums are needed for step eq. (forward and backward) and for cross-complementarties)
            Lambda_ki_all = [Lambda_end_previous_fe,Lambda_ki];
        end
        % Update the standard complementarity
        % TODO!!: how is this function evaluated properly with a matrix?
        J_comp_std = J_comp_std + J_cc_fun(Z_ki);
        %% Additional boundary points if c_{n_s} \neq 1. (cf. FESD paper)
%         if use_fesd && ~right_boundary_point_explicit &&  (k<N_stages-1 || i< N_finite_elements(k+1)-1)
%             switch pss_mode
%                 case 'Stewart'
%                     %                     eval(['Lambda_ki_end  = ' casadi_symbolic_mode '.sym(''Lambda_' num2str(k) '_' num2str(i) '_end' ''',n_theta);'])
%                     %                     eval(['Mu_ki_end  = ' casadi_symbolic_mode '.sym(''Mu_' num2str(k) '_' num2str(i) '_end' ''',n_simplex);'])
%                     Lambda_ki_end = define_casadi_symbolic(casadi_symbolic_mode,['Lambda_' num2str(k) '_' num2str(i) '_end'],n_theta);
%                     Mu_ki_end = define_casadi_symbolic(casadi_symbolic_mode,['Mu_' num2str(k) '_' num2str(i) '_end'],n_simplex);
%                     w = {w{:}, Lambda_ki_end};
%                     w = {w{:}, Mu_ki_end};
%                     % bounds and index sets
%                     w0 = [w0; z0(n_theta+1:end)];
%                     lbw = [lbw; lbz(n_theta+1:end)];
%                     ubw = [ubw; ubz(n_theta+1:end)];
%                     ind_boundary = [ind_boundary,ind_total(end)+1:(ind_total(end)+n_z-n_theta)];
%                     ind_total  = [ind_total,ind_total(end)+1:(ind_total(end)+n_z-n_theta)];
%                     Lambda_ki_current_fe = {Lambda_ki_current_fe{:}, Lambda_ki_end};
%                     Mu_ki_current_fe = {Mu_ki_current_fe{:}, Mu_ki_end};
%                 case 'Step'
%                     %                     eval(['Lambda_ki_end  = ' casadi_symbolic_mode '.sym(''Lambda_' num2str(k) '_' num2str(i) '_end' ''',2*n_alpha);'])
%                     Lambda_ki_end = define_casadi_symbolic(casadi_symbolic_mode,['Lambda_' num2str(k) '_' num2str(i) '_end'],2*n_alpha);
%                     Mu_ki_end = [];
%                     w = {w{:}, Lambda_ki_end};
%                     % bounds and index sets
%                     w0 = [w0; z0(n_alpha+1:3*n_alpha)];
%                     lbw = [lbw; lbz(n_alpha+1:3*n_alpha)];
%                     ubw = [ubw; ubz(n_alpha+1:3*n_alpha)];
%                     ind_boundary = [ind_boundary,ind_total(end)+1:(ind_total(end)+2*n_alpha)];
%                     ind_total  = [ind_total,ind_total(end)+1:(ind_total(end)+2*n_alpha)];
%                     Lambda_ki_current_fe = {Lambda_ki_current_fe{:}, Lambda_ki_end};
%                     Mu_ki_current_fe = {Mu_ki_current_fe{:}, []};
%             end
%             % For cross comp
%             Lambda_sum_finite_element_ki =  Lambda_sum_finite_element_ki + Lambda_ki_end;
%         end

        %% Continuity of Lambda - Reinitalize Lambda_end_previous_fe
%             !! TODO
 
        % Defintion of all stage variables for current finite element is
        % done]
        %% The IRK Equations - evaluation of dynamics (ODE r.h.s. and algebraic equations at stage points)
        
       
       switch irk_representation
            case 'integral'
                % Initalize sum for quadrature of 
%                 Xk_end = D(1)*X_ki; <--  TODO: Note this expressions is not used anymore
                
                % Note that the polynomial is initalized with the previous value % (continuity connection)
                % Xk_end = D(1)*Xk;   % X_k+1 = X_k0 + sum_{j=1}^{d} D(j)*X_kj; Xk0 = D(1)*Xk
                % X_k+1 = X_k0 + sum_{j=1}^{n_s} D(j)*X_kj; Xk0 = D(1)*Xk
                if n_u > 0
                    [f_x, f_q] = f_x_fun(X_ki_stages, Z_ki, U_k);
                    [g_z] = g_lp_fun(X_ki_stages, Z_ki, U_k);
                else
                    [f_x, f_q] = f_x_fun(X_ki_stages, Z_ki);
                    [g_z] = g_lp_fun(X_ki_stages, Z_ki);
                    % TODO: distignush genral g_lp_fun and g_z
                end

                % Time-rescaling of dynamics
                if time_rescaling && use_speed_of_time_variables
%                         % rescale equations
%                         fj = s_sot_k*fj;
%                         qj = s_sot_k*qj;
%                         gj = s_sot_k*gj;
                end
                % All stage values in one matrix 
                X_all = [X_ki X_ki_stages];
                if use_fesd
                    J = J + f_q*B*h_ki;
                    % Get interpolating points of collocation polynomial
                    
                    % Get slope of interpolating polynomial (normalized)
                    Pidot = X_all*C;
                    % Match with ODE right-hand-side
                    opti.subject_to(Pidot == h_ki*f_x);
                    opti.subject_to(0 == g_z);
                else
                    J = J + f_q*B*h;
                    % Get interpolating points of collocation polynomial
                    X_all = [X_ki X_ki_stages];
                    % Get slope of interpolating polynomial (normalized)
                    Pidot = X_all*C;
                    % Match with ODE right-hand-side
                    opti.subject_to(Pidot == h*f_x);
                    opti.subject_to(0 == g_z);
                end
                % State at end of finite elemnt
                Xk_end = X_all*D;

            case 'differential'
                Xk_end = X_ki; % initalize with x_n;
       end

       %% Box constraint at stage points if now lifting is done (via expressions for X_ki_j) in differnatil mode
        % TODO ; do this above where they are defined

       %% Lifting expressions
        % Todo : be careful when using a matrix expression in subject_to
    %                if lift_irk_differential
    %                 g = {g{:}, X_ki_lift{j}};
    %                 lbg = [lbg; zeros(n_x,1)];
    %                 ubg = [ubg; zeros(n_x,1)];
    %             end
       %% General nonlinear constraint at stage point 
        if g_ineq_constraint && g_ineq_at_stg

        end

        %% Cross complementarity and step equlibration constraints

        %% Treatment of complementarity constraints (MPCC reformulation)

        %% Continuity condition for differential state - new NLP variable for state at end of a finite element
        % New decision variable for state at end of a finite element
        X_ki = opti.variable(n_x);
        X_boundary{end+1} = X_ki;
        X{end+1} = X_ki;
        opti.subject_to(lbx <= X_ki <= ubx);
        opti.set_initial(X_ki, x0);
        ind_x= [ind_x,ind_total(end)+1:ind_total(end)+n_x];
        ind_total = [ind_total,ind_total(end)+1:ind_total(end)+n_x];
        % Continuity constraints
        opti.subject_to(Xk_end==X_ki)
        %% Evaluate inequiality at finite element boundary?
         % TODO how does this compile with the evaulation after defining the control variables

        %% G_LP constraint for boundary point and continuity of algebraic variables.
          % TODO: give a cleaner discription, this the partial algebraic eq. for c_{n_s} \neq 1 
    end
    sum_h_ki = [sum_h_ki;sum_h_ki_temp];
    %% Equdistant grid in numerical time (Multiple-shooting type discretization)

    %% Equdistant grid in phyisical time (Stage-wise constraints on the colock state)

end
sum_h_ki_all = sum(sum_h_ki); % this is the sum of all over all FE and control intervals (= integral of clock state if no time-freezing is used)
%% Terminal Constraints
%  possible relaxation of terminal constraints

% terminal constraint for physical and numerical  time

%% Terminal Costs
% basic terminal cost term

% Quadrature states

% Time-optimal problems

%% Elastic mode cost terms

%% Regularization cost terms
% Regularization term for speed-of-time
% TODO: Idea: add sot to model reformulation and add f_q (or do it all here, to
% Cost term for grid regularization (step equilibration, heursitic step)

%% Collect all variables
X = [X{:}];
X_boundary = [X_boundary{:}];
V = [V{:}];
U = [U{:}];
Z = [Z{:}];
H = [H{:}];
S_sot = [S_sot{:}];
S_elastic = [S_elastic{:}];
%% CasADi functions for complementarity residuals (standard, cross_complementarity, joint)

%% Create NLP Solver instance
opti.minimize(J);
opti.solver('ipopt'); % TODO: IPOPT Settings;
% opti.solver(solver_name,'ipopt',opts_ipopt);
%TODO! : Loop of problems for incerasing number of iterations;
% e.g. min iter = 100, max_iter, 1500, make ~equidistant grid of integers on ceil(linspace(min_iter,max_iter,N_homotopy); save all solvers in cell
% sol = opti.solve();
%% Results
% x_opt = sol.value(X_boundary);
% x_opt_extended = sol.value(X);
% z_opt = sol.value(Z);
% u_opt = sol.value(U);
% J_opt = sol.value(J);

%% Model: CasADi functions, dimenesions, auxilairy functions.

%% Solve initalization (bounds, inital guess, parameters)

%% Output of the function
varargout{1} = solver;
varargout{2} = solver_initalization;
varargout{3} = model;
varargout{4} = settings;
end