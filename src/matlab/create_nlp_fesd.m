function [solver,solver_initalization, model,settings] = create_nlp_fesd(model,settings)
% This functions creates the solver instance for the OCP discretized with FESD (or time stepping IRK scheme).
% The discretization results in an MPCC which can be solved by various
% reformulations, see below.
% -------------------
% Brief MPCC Wiki
% There are several possible MPCC Solution strategies avilable, by setting mpcc_mode to :
% 1 - treat complementarity conditions directly in the NLP, the bilinear term is tread as an inequality constraint.
% 2 - Smooth the complementarity conditions.
% 3 - Relax the complementarity conditions.
% 4 - \ell_1 penalty, penalize the sum of all bilinear terms in the objective
% 5 - \ell__infty elastic mode, upper bound all bilinear term with a positive slack, and penalize the slack in the objective.
% 6 - \ell__infty elastic mode, equate all bilinear term to a positive slack, and penalize the slack in the objective.
% 7 - \ell__infty, same as 5 but two sided.
% 8 - \ell__infty, same as 5 but the penalty/slack is controlled via the barier formulation
% 9 - \ell__infty, same as 6 but the penalty/slack is controlled via the barier formulation
% 10 - \ell__1, same as 4 but the penalty/slack is controlled via the barier formulation
% Main developer: Armin Nurkanovi\'c (armin.nurkanovic@imtek.uni-freiburg.de)
%% Import CasADi in the workspace of this function
import casadi.*

%% Reformulation of the PSS into a DCS
[model,settings] = model_reformulation_fesd(model,settings);
%% Fillin missing settings with default settings
[settings] = fill_in_missing_settings(settings,model);

%% Load user settings and model details
unfold_struct(settings,'caller')
% Note that some default values of N_stages ect are overwritten by the model data, if provided.
unfold_struct(model,'caller');

%% Little corrections
% TODO: where should i put the model refinment code ...
if length(N_finite_elements) == 1
    N_finite_elements = N_finite_elements*ones(N_stages,1);
elseif length(N_finite_elements) > 1 && length(N_finite_elements) < N_stages
    N_finite_elements = N_finite_elements(:); % make sure it is a column vector
    N_finite_elements = [N_finite_elements;N_finite_elements(end)*ones(N_stages-length(N_finite_elements),1)];
end

h_k = h./N_finite_elements; %(TODO:take care if this is a vector)
%% Create a structure with dimensions, useful for some other functions
%% TODO move this in one of the refinemnt functions
dimensions.N_stages = N_stages;
dimensions.N_finite_elements = N_finite_elements;
dimensions.n_x = n_x;
dimensions.n_u = n_u;
dimensions.n_z = n_z;
dimensions.d= d;
dimensions.n_simplex = n_simplex;
dimensions.m_vec = m_vec;
dimensions.m_ind_vec = m_ind_vec;

%% Elastic Mode Variables
if mpcc_mode >= 5 && mpcc_mode <= 10
    s_elastic = MX.sym('s_elastic', 1);
end
%% Initalization and bounds
% solve LP for guess;
if lp_initalization
    [theta_guess,lambda_guess,mu_guess] = create_lp_based_guess(model);
else
    theta_guess = initial_theta*ones(n_theta,1);
    lambda_guess = initial_lambda*ones(n_theta,1);
    mu_guess = initial_mu*ones(n_simplex,1);
end
z0 = [theta_guess;lambda_guess;mu_guess];
% Lower and upper bounds for \theta, \lambda and \mu.
lbz = [0*ones(n_theta,1);0*ones(n_theta,1);-inf*ones(n_simplex,1)];
ubz = [inf*ones(n_theta,1);inf*ones(n_theta,1);inf*ones(n_simplex,1)];

%% Initalization and bounds of algebraic variables
if use_fesd
    lbh = (1-gamma_h)*h_k;
    ubh = (1+gamma_h)*h_k;

    if time_rescaling && ~use_speed_of_time_variables
        % if only time_rescaling is true, speed of time and
        % step size all lumped together, e.g., \hat{h}_{n,m} = s_n * h_{n,m}, hence the bounds need to be extended.
        ubh = (1+gamma_h)*h_k*s_sot_max;
        lbh = (1-gamma_h)*h_k/s_sot_min;
    end
    % initigal guess for the step-size
    h0_k = h_k;
end

%% Time optimal control problems without speed of time variables
% if time_optimal_problem && (~use_speed_of_time_variables || time_freezing)
if time_optimal_problem
    % this should be independet weathe or not time freezing is used.
    T_final = MX.sym('T_final',1);
    T_final_guess = T;
end

%% Collocation related values, Butcher Tableu
[B,C,D,tau_root] = collocation_times_fesd(d,collocation_scheme);
%% Formulate NLP
% Start with an empty NLP
w = {};
w0 = [];
lbw = [];
ubw = [];
J = 0;
J_comp = 0;
J_regularize_h = 0;
J_regularize_h_abs = 0;
J_regularaize_delta_h  = 0;
g={};
lbg = [];
ubg = [];

% group the constraints as in Section 5 of the FESD Paper.
g_irk = {};
lbg_irk = [];
ubg_irk = [];

g_cross = {};
lbg_cross = [];
ubg_cross = [];

g_cont = {};
lbg_cont = [];
ubg_cont = [];

g_eq = {};
lbg_eq = [];
ubg_eq = [];

% "Lift" initial conditions
X_kq = MX.sym('X0', n_x);
w = {w{:}, X_kq};

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
ind_z = [];
ind_h = [];
ind_g_clock_state = [];
ind_tf = []; % time freezing index
ind_boundary = []; % index of bundary value lambda and mu
ind_total = ind_x;

sum_h_kq_all_stages = 0;
sum_h_kq_stage_k = 0;

% Initalization of forward and backward sums for the step-equilibration

sigma_lambda_F_k = 0; % forward sum of lambda at stage k
sigma_lambda_B_k = 0; % backward sum of lambda at stage k

sigma_theta_F_k = 0; % forward sum of theta at stage k
sigma_theta_B_k = 0; % backward sum of theta at stage k

nu_vector = [];

% Integral of the clock state if no time freezing is present.
sum_of_S_sot_k =0;
n_cross_comp = zeros(max(N_finite_elements),N_stages);

%% Formulate the NLP / Main Discretization loop
for k=0:N_stages-1
    %% NLP variables for the controls
    Uk = MX.sym(['U_' num2str(k)],n_u);
    w = {w{:}, Uk};
    % intialize contros, lower and upper bounds
    w0 = [w0; u0];
    lbw = [lbw;lbu];
    ubw = [ubw;ubu];
    % index colector for contorl variables
    ind_u = [ind_u,ind_total(end)+1:ind_total(end)+n_u];
    ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_u];

    %%  Time rescaling of the stages (speed of time) to acchieve e.g., a desired final time in Time-Freezing or to solve time optimal control problems.
    %     If only time_rescaling is true, then the h_k also make sure to addapt the length of the subintervals, if both
    %     time_rescaling && use_speed_of_time_variables are true, new control variables are introduecd, they can be per stage or one for the whole interval.
    if time_rescaling && use_speed_of_time_variables
        if local_speed_of_time_variable
            % at every stage
            S_sot_k = MX.sym(['S_sot_' num2str(k)],1);
            w = {w{:}, S_sot_k};
            % intialize speed of time (sot), lower and upper bounds
            w0 = [w0; s_sot0];
            ubw = [ubw;s_sot_max];
            lbw = [lbw;s_sot_min];
            % index colector for sot variables
            ind_tf = [ind_tf,ind_total(end)+1:ind_total(end)+1];
            ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];
        else
            % only once
            if k == 0
                S_sot_k = MX.sym(['S_sot_' num2str(k+1)],1);
                w = {w{:}, S_sot_k};
                % intialize speed of time (sot), lower and upper bounds
                w0 = [w0; s_sot0];
                lbw = [lbw;s_sot_min];
                ubw = [ubw;s_sot_max];
                % index colector for sot variables
                ind_tf = [ind_tf,ind_total(end)+1:ind_total(end)+1];
                ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];
            end
        end
    end

    %% General Nonlinear constriant (on finite element boundary)
    % The CasADi function g_ineq_fun and its lower and upper bound are provieded in model.
    if general_nonlinear_constraint
        g_ineq_k = g_ineq_fun(X_kq,Uk(1:n_u));
        g = {g{:}, g_ineq_k };
        lbg = [lbg; g_ineq_lb];
        ubg = [ubg; g_ineq_ub];
    end
    sum_h_kq_stage_k = 0; % Integral of the clock state (if not time freezing) on the current stage.

    %% Define step-size variables, if FESD is used.
    for i = 0:N_finite_elements(k+1)-1
        if use_fesd
            if  k>0 || i>0
                h_kq_previous = h_kq;
            end
            % define step-size
            h_kq = MX.sym(['h_' num2str(k) '_' num2str(i)],1);

            w = {w{:}, h_kq};
            w0 = [w0; h0_k(i+1)];
            ubw = [ubw;ubh(i+1)];
            lbw = [lbw;lbh(i+1)];

            % index colector for contorl variables
            ind_h = [ind_h,ind_total(end)+1:ind_total(end)+1];
            ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];

            sum_h_kq_stage_k = sum_h_kq_stage_k + h_kq;
            sum_h_kq_all_stages = sum_h_kq_all_stages + h_kq;

            if time_optimal_problem && use_speed_of_time_variables
                sum_of_S_sot_k = sum_of_S_sot_k + h_kq*S_sot_k;
            end

            if (k > 0 && (i > 0 || couple_across_stages))
                delta_h_kq = h_kq - h_kq_previous;
            else
                delta_h_kq  = 0;
            end

            if heuristic_step_equilibration && (k > 0 && (i > 0 || couple_across_stages))
                % quadratic term for getting equidistant grid
                switch heuristic_step_equilibration_mode
                    case 1
                        J_regularize_h  = J_regularize_h + (h_kq-h_k(i+1))^2;

                    case 2
                        J_regularize_h  = J_regularize_h + delta_h_kq^2;
                    otherwise
                        error('Pick heuristic_step_equlibration_mode between 1 and 2.');
                end
            end

            %% The crosscomplementartiy condtions: untilize sums of \lambda and \theta
            % \sum_{j=0}^{d} \lambda_{k,j}  amd % \sum_{i=1}^{d} \theta_{k,i}
            % calculate algerbaric variale at t = t0,(needed for complementarity quadrature;)
            %% TODO: Does this already inlcude coupling between all finite elements??
            if k == 0 && i == 0
                Z_kdq = zeros(n_z,1);
            end
            % initalize with last collocation point from previous finite element
            sum_lambda_ki = Z_kdq(n_theta+1:2*n_theta);
            % initalize sum of theta's (the pint at t_k is not included)
            sum_theta_ki = 0;
        end

        %% Define Variables at Collocation Points (IRK stages) at all Finite Elements
        % State at collocation points
        X_kqj = {};
        Z_kqj = {}; % collects the vector of x and z at every collocation point/irk stage
        if i>0 && use_fesd 
            Lambda_end_previous_fe = {Z_kdq(n_theta+1:2*n_theta)};
        end
        Lambda_ki_this_fe = {};
        Theta_ki_this_fe = {};
        Mu_ki_this_fe = {};

        for j=1:d
            X_kqj{j} = MX.sym(['X_' num2str(k) '_' num2str(i) '_' num2str(j)], n_x);
            Z_kqj{j} = MX.sym(['Z_' num2str(k) '_' num2str(i) '_' num2str(j)], n_z);

            w = {w{:}, X_kqj{j}};
            w = {w{:}, Z_kqj{j}};

            % Index collector for algebraic and differential variables
            ind_x = [ind_x,ind_total(end)+1:ind_total(end)+n_x];
            ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_x];

            ind_z = [ind_z,ind_total(end)+1:ind_total(end)+n_z];
            ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_z];

            % Bounds
            lbw = [lbw; lbx];
            ubw = [ubw; ubx];

            lbw = [lbw; lbz];
            ubw = [ubw; ubz];

            % Trajectory initial guess
            w0 = [w0; x0];
            w0 = [w0; z0];

            if use_fesd
                % sum_{j = 0}^{d} lambda_{k,j} (sum lambdas over all collocation points)
                sum_lambda_ki =  sum_lambda_ki + Z_kqj{j}(n_theta+1:2*n_lambda);
                % sum_{i = 1}^{d} lambda_{k,i} (sum thetas over all collocation points)
                sum_theta_ki =  sum_theta_ki +  Z_kqj{j}(1:n_theta);
            end

            % collection of all lambda and theta for current finite element (stage), \lambda_{n,q}, \theta_{n,q} in the paper in Section 4 (or \lambda_{n}  and \theta_{n} in other sections

            Theta_ki_this_fe = {Theta_ki_this_fe{:}, Z_kqj{j}(1:n_theta)};
            Lambda_ki_this_fe = {Lambda_ki_this_fe{:}, Z_kqj{j}(n_theta+1:2*n_lambda)};
            Mu_ki_this_fe = {Mu_ki_this_fe{:}, Z_kqj{j}(end)};
        end
        n_cross_comp_i = 0;
        g_cross_comp_i = []; % all cross complementarites over a finite element

        %% Additiona variables in case of schemes not containting the boundary point as collocation point, e.g., Gauss-Legendre schemes
        g_cont_i = [];
        if use_fesd 
%             &&  isequal(collocation_scheme,'legendre')
            % define variables for terminal values
            Lambda_ki_end = MX.sym(['Lambda_' num2str(k) '_' num2str(i) '_end'], n_theta);
            Mu_ki_end = MX.sym(['Mu_' num2str(k) '_' num2str(i) '_end'], 1);

            w = {w{:}, Lambda_ki_end};
            w = {w{:}, Mu_ki_end};

            Mu_ki_this_fe_vec = [];
            Lambda_ki_this_fe_vec = [];
            for j = 1:d
                Mu_ki_this_fe_vec = [Mu_ki_this_fe_vec,Mu_ki_this_fe{j}];
                Lambda_ki_this_fe_vec = [Lambda_ki_this_fe_vec,Lambda_ki_this_fe{j}];
            end

            g_cont_i = [g_cont_i;Mu_ki_end-Mu_ki_this_fe_vec*m_lagrange'];
            g_cont_i = [g_cont_i;Lambda_ki_end-Lambda_ki_this_fe_vec*m_lagrange'];


            g = {g{:}, g_cont_i};
            lbg = [lbg; zeros(n_theta+1,1)];
            ubg = [ubg; zeros(n_theta+1,1)];

            %                 % Index collector for algebraic and differential variables
            ind_boundary = [ind_boundary,ind_total(end)+1:(ind_total(end)+n_theta+1)];
            ind_total  = [ind_total,ind_total(end)+1:(ind_total(end)+n_theta+1)];

            %                 % Bounds
%             lbw = [lbw; -inf*ones(n_theta+1,1)];
            lbw = [lbw; [0*ones(n_theta,1);-inf*ones(1,1)]];
            ubw = [ubw; +inf*ones(n_theta+1,1)];
            %                % Trajectory initial guess
            w0 = [w0; z0(n_theta+1:end)];
            if isequal(collocation_scheme,'legendre')
                sum_lambda_ki =  sum_lambda_ki + Lambda_ki_end;
            end
        end
        
        %% Evaluate equations (dynamics, algebraic, complementarities standard and cross at every collocation point)
        % Note that the polynomial is initalized with the previous value % (continuity connection)
        % Xk_end = D(1)*Xk;   % X_k+1 = X_k0 + sum_{j=1}^{d} D(j)*X_kj; Xk0 = D(1)*Xk
        % X_k+1 = X_k0 + sum_{j=1}^{d} D(j)*X_kj; Xk0 = D(1)*Xk
        Xk_end = D(1)*X_kq;
        g_cross_comp_temp = 0; % temporar sum for cross comp (whenever there is a single constraint per finite element, casese 4 and 6)
        for j=1:d
            % ODE r.h.s. and 'standard' algebraic equations
            % Expression for the state derivative at the collocation point
            xp = C(1,j+1)*X_kq;

            % Lagrange polynomial with values at state
            for r=1:d
                xp = xp + C(r+1,j+1)*X_kqj{r};
            end
            % Evaulate Differetinal and Algebraic Equations at Collocation Points
            [fj, qj] = f_x(X_kqj{j},Z_kqj{j},Uk(1:n_u));
            gj = f_z(X_kqj{j},Z_kqj{j},Uk(1:n_u));

            if time_rescaling && use_speed_of_time_variables
                % rescale equations
                fj = S_sot_k*fj;
                qj = S_sot_k*qj;
                gj = S_sot_k*gj;
            end

            % Append collocation equations. append differential states equations and algebraic complementarity conditions
            if use_fesd
                g = {g{:}, h_kq*fj - xp};
            else
                g = {g{:}, h_k(k+1)*fj - xp};
            end
            g = {g{:}, gj};
            lbg = [lbg; zeros(n_x,1); zeros(n_algebraic_constraints,1)];
            ubg = [ubg; zeros(n_x,1); zeros(n_algebraic_constraints,1)];

            % Add contribution to the end state
            Xk_end = Xk_end + D(j+1)*X_kqj{j};
            % Add contribution to quadrature function
            if use_fesd
                J = J + B(j+1)*qj*h_k(i+1);
            else
                J = J + B(j+1)*qj*h;
            end

            % Collcet all IRK equations:
            if use_fesd
                g_irk = {g_irk{:} , h_kq*fj - xp};
            else
                g_irk = {g_irk{:} , h_k(k+1)*fj - xp};
            end
            g_irk = {g_irk{:}, gj};
            lbg_irk = [lbg_irk; zeros(n_x,1); zeros(n_algebraic_constraints,1)];
            ubg_irk = [ubg_irk; zeros(n_x,1); zeros(n_algebraic_constraints,1)];


            %% Standard and cross complementarity constraints
            % Prepare Input for Cross Comp Function
            comp_var_this_fe.J_comp = J_comp;
            try
                comp_var_this_fe.sum_lambda_ki = sum_lambda_ki;
                comp_var_this_fe.sum_theta_ki = sum_theta_ki;
            catch
                % not avilable if not use_fesd = 0. TODO: maybe still have
                % it for different sparsity in the use_fesd =0 case.
            end
            comp_var_this_fe.Lambda_ki_this_fe = Lambda_ki_this_fe;
            comp_var_this_fe.Theta_ki_this_fe = Theta_ki_this_fe;
            current_index.k = k;  current_index.i = i; current_index.j = j;

            % Create cross comp constraints
            [J_comp,g_cross_comp_j] = create_complementarity_constraints(use_fesd,cross_complementarity_mode,comp_var_this_fe,dimensions,current_index);

            n_cross_comp_j = length(g_cross_comp_j);
            n_cross_comp_i = n_cross_comp_i+n_cross_comp_j;
%             g_cross_comp_i = [g_cross_comp_i;g_cross_comp_j]; % IS THIS even used anywhere? Should be used for storing of FESD equations?

            n_cross_comp(i+1,k+1) = n_cross_comp_i;
            g_all_comp_j = [g_cross_comp_j];
            n_all_comp_j = length(g_all_comp_j);

            %% Treatment and reformulation of all Complementarity Constraints (standard and cross complementarity), their treatment depends on the chosen MPCC Method.
            mpcc_var_this_fe.J = J;
            mpcc_var_this_fe.g_all_comp_j = g_all_comp_j;
            % TODO: this could be ented just once at the beggining as they dont change
            mpcc_var_this_fe.p = p;
            try
                mpcc_var_this_fe.s_elastic = s_elastic;
            catch
                % this makes sense only if s_elastic exists
            end
            [J,g_comp,g_comp_lb,g_comp_ub] = reformulate_mpcc_constraints(objective_scaling_direct,mpcc_mode,mpcc_var_this_fe,dimensions,current_index);
            g = {g{:}, g_comp};
            lbg = [lbg; g_comp_lb];
            ubg = [ubg; g_comp_ub];

            %% General nonlinear constraint at colocation points
            if general_nonlinear_constraint && general_nonlinear_constraint_at_collocation_points
                g_ineq_k = g_ineq_fun(X_kqj{j},Uk(1:n_u));
                g = {g{:}, g_ineq_k};
                lbg = [lbg; g_ineq_lb];
                ubg = [ubg; g_ineq_ub];
            end
        end

        %% Continuity conditions
        if use_fesd
            switch collocation_scheme
                case 'radau'
                    Z_kdq = Z_kqj{d};
                case 'lobbato'
                      Z_kdq = Z_kqj{d};
                case 'legendre'
                    Z_kdq = [zeros(n_theta,1);Lambda_ki_end;Mu_ki_end];
                otherwise
                    error('Pick Gauss-Legendre or Radau schemes.');
            end
        end

        %% Step equlibration
        % Define the functions for propper grid regularization, sigma, delta, pi ...
        % TODO: MAKE A Function of the script below
        step_equilibration_constrains;

        %% Continuity condition  New NLP variable for state at end of a finite element
        % TODO: This the novel variable for the contuinity, should be renamed such that it maches the other notation, cf defintion of collocation points X_kq --> X_kq0
        X_kq = MX.sym(['X_' num2str(k+1) '_' num2str(i+1)], n_x);
        w = {w{:}, X_kq};

        ind_x = [ind_x,ind_total(end)+1:ind_total(end)+n_x];
        ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+n_x];

        lbw = [lbw; lbx];
        ubw = [ubw; ubx];
        w0 = [w0; x0];

        % Add equality constraint
        g = {g{:}, Xk_end-X_kq};
        lbg = [lbg; zeros(n_x,1)];
        ubg = [ubg; zeros(n_x,1)];

        %% Switch detecting constraint
        if isequal(collocation_scheme,'legendre') && k< N_stages-1
            gj = f_z(X_kq,Z_kdq,Uk(1:n_u));
            g = {g{:}, gj(1:end-1)};
            lbg = [lbg; zeros(n_theta,1)];
            ubg = [ubg; zeros(n_theta,1)];
        end
    end

    %% Equdistant grid in numerical time (Multiple-shooting type discretization)
    if equidistant_control_grid
        if time_rescaling && time_optimal_problem
            % the querry here is because: No Time freezing: time_opt => time_rescaling (so this is always true if time_rescaling is on)
            % If time freezing is present, but not ime optimal problem, then the final numerical time should not be changed, hence the query.
            if use_speed_of_time_variables
                % numerical time and speed of time are decoupled
                g = {g{:}, sum_h_kq_stage_k-h};
                lbg = [lbg; 0];
                ubg = [ubg; 0];
            else
                % numerical time and speed of time are lumped together
                g = {g{:}, sum_h_kq_stage_k-T_final/N_stages};
                lbg = [lbg; 0];
                ubg = [ubg; 0];
            end
        else
            % numerical time is fixed.
            g = {g{:}, sum_h_kq_stage_k-h};
            lbg = [lbg; 0];
            ubg = [ubg; 0];
        end
    end

    %% Equdistant grid in phyisical time (Stage-wise constraints on the colock state)
    if time_freezing && stagewise_clock_constraint
        % This makes mostly sense if time freezin is on. Imposes the constraints t((k+1)*h) = (k+1)*h+t_0 , x0(end) = t_0.
        if time_optimal_problem
            % this seems to be a critical constraint, T_final has to be
            % set equal to the sum_sot at the end?
            g = {g{:}, Xk_end(end)-(k+1)*(T_final/N_stages)+x0(end)};
        else
            g = {g{:}, Xk_end(end)-(k+1)*h+x0(end)};
        end
        lbg = [lbg; 0];
        ubg = [ubg; 0];
        ind_g_clock_state = [ind_g_clock_state;length(lbg)];
    end
end

%% Scalar-valued commplementarity residual
J_comp = sum(J_comp);
%% Constraint for the terminal numerical and physical time (if no equidistant grids are required)
% If the control grid is not equidistant, the constraint on sum of h happen only at the end.
% The constraints are splited to those which are related to numerical and physical time, to make it easier to read.

%% Numerical Time
if  ~equidistant_control_grid
    if time_rescaling
        % this means this is a time optimal control problem without time freezing
        if use_speed_of_time_variables
            % numerical time is decupled from the sot scaling (no mather if local or not):
            g = {g{:}, sum_h_kq_all_stages-T};
            lbg = [lbg; 0];
            ubg = [ubg; 0];
            % T = T_num,  T_phy = s_sot*T.
            % note that objective contribution is added below (as it happens no
            % mather if we use local or nonlocal variables and equdistiant grid or not.)
        else
            if time_optimal_problem
                % the querry here is because: No Time freezing: time_opt => time_rescaling (so this is always true if time_rescaling is on)
                % If time freezing is present, but not ime optimal problem, then the final numerical time should not be changed, hence the query.
                g = {g{:}, sum_h_kq_all_stages-T_final};
                lbg = [lbg; 0];
                ubg = [ubg; 0];
                % add time to objective.
                J = J+T_final;
            end
            % T_num = T_final = T_phy \neq T.
        end
    else
        % fixed numerical time (no time freezing or time optimal control, T_phy = T_num = T)
        g = {g{:}, sum_h_kq_all_stages-T};
        lbg = [lbg; 0];
        ubg = [ubg; 0];
    end
end

%% Phyisical Time (Posssble terminal constraint on the clock state if time freezing is active).
if time_freezing
    if time_optimal_problem
        % NOT TESTED
        g = {g{:}, Xk_end(end)-T_final};
        lbg = [lbg; 0];
        ubg = [ubg; 0];
    else
        if impose_terminal_phyisical_time && ~stagewise_clock_constraint
            % THOS
            g = {g{:}, Xk_end(end)-T};
            lbg = [lbg; 0];
            ubg = [ubg; 0];
        else
            % no terminal constraint on the numerical time, as it
            % is implicityl determined by
        end
    end
end
% TODO: FInd out what is this covering?
if equidistant_control_grid && ~stagewise_clock_constraint && time_freezing
    if ~time_optimal_problem
        g = {g{:}, Xk_end(end)-T};
        lbg = [lbg; 0];
        ubg = [ubg; 0];
    else

    end
end

%%  Time-optimal problems: the objective.

if time_optimal_problem && use_speed_of_time_variables
    if time_freezing
        J = J+T_final;
    else
        % This is both variants local and global sot vars. Hence, it is written here together and not in the loops or queries above.
        % ( In this case,the trivial constraint T_final = sum_of_S_sot_k) is ommited and provided implicitly)
        g = {g{:}, sum_of_S_sot_k-T_final};
        lbg = [lbg; 0];
        ubg = [ubg; 0];
        J = J+T_final;
        %         J = J+sum_of_S_sot_k;
    end
end

% Time-optimal problems: define auxilairy variable for the final time.
ind_t_final = [];
% if time_optimal_problem && (~use_speed_of_time_variables || time_freezing)
if time_optimal_problem
    % Add to the vector of unknowns
    w = {w{:}, T_final};
    w0 = [w0; T_final_guess];
    lbw = [lbw;0];
    ubw = [ubw;1e2];
    ind_t_final = [ind_total(end)+1];
    ind_total  = [ind_total,ind_total(end)+1];
end

%% Terminal Constraints
% Add Terminal Constrint
if terminal_constraint
    g_terminal = g_terminal_fun(Xk_end);
    g = {g{:}, g_terminal };
    lbg = [lbg; g_terminal_lb];
    ubg = [ubg; g_terminal_ub];
end

%% Add Terminal Cost to objective
try
    J = J + f_q_T(Xk_end);
catch
    warning('Terminal cost not defined');
end
J_objective = J;

%% Elastic mode variable for \ell_infty reformulations
if mpcc_mode >= 5 && mpcc_mode < 8
    % add elastic variable to the vector of unknowns and add objeective contribution
    w = {w{:}, s_elastic};
    ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];
    lbw = [lbw; s_elastic_min];
    ubw = [ubw; s_elastic_max];
    w0 = [w0;s_elastic_0];

    if objective_scaling_direct
        J = J + (1/p)*s_elastic;
    else
        J = p*J + s_elastic;
    end
end

%% barrier controled penalty formulation
if mpcc_mode >= 8 && mpcc_mode <= 10
    rho_elastic = MX.sym('rho_elastic', 1);
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

    % add elastic variable to the vector of unknowns and add objeective contribution
    w = {w{:}, s_elastic};
    ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];
    lbw = [lbw; s_elastic_min];
    ubw = [ubw; s_elastic_max];
    w0 = [w0;s_elastic_0];

    J = J-rho_penalty*(rho_elastic^2)+sigma_penalty*s_elastic;
end

%% Objective Terms for Grid Regularization
% Huristic Regularization.
if heuristic_step_equilibration
    J = J + step_equilibration_penalty*J_regularize_h;
end

if step_equilibration
    J = J + step_equilibration_penalty*J_regularize_h;
end
%% CasADi Functions for objective complementarity residual
J_fun = Function('J_fun', {vertcat(w{:})},{J_objective});
comp_res = Function('comp_res',{vertcat(w{:})},{J_comp});

%% NLP Solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}),'p',p);
solver = nlpsol(solver_name, 'ipopt', prob,opts_ipopt);

%% Define CasADi function for the switch indicator function.
if step_equilibration
    nu_fun = Function('nu_fun', {vertcat(w{:})},{nu_vector});
end

%% Outputs
model.prob = prob;
model.solver = solver;
model.g =  vertcat(g{:});
model.w =  vertcat(w{:});
model.J = J;
model.J_fun = J_fun;
model.comp_res = comp_res;

% create CasADi function for objective gradient.
nabla_J = J.jacobian(model.w);
nabla_J_fun = Function('nabla_J_fun', {vertcat(w{:}),p},{nabla_J});

model.nabla_J = nabla_J;
model.nabla_J_fun = nabla_J_fun;

if step_equilibration
    model.nu_fun = nu_fun;
end

%% Store solver initalization data
solver_initalization.w0 = w0;
solver_initalization.lbw = lbw;
solver_initalization.ubw = ubw;
solver_initalization.lbg = lbg;
solver_initalization.ubg = ubg;
n_cross_comp_total = sum(n_cross_comp(:));

%% Model update: all index sets and dimensions
model.ind_x = ind_x;
model.ind_z = ind_z;
model.ind_u = ind_u;
model.ind_h = ind_h;
model.ind_g_clock_state = ind_g_clock_state;
model.ind_tf = ind_tf;
model.ind_t_final  = ind_t_final;
model.ind_boundary = ind_boundary;

model.n_u = n_u; % dimenison of control without variable step size
model.n_cross_comp = n_cross_comp;

model.ind_total = ind_total;
model.nlp = prob;

settings.h = h;
settings.h_k = h_k;
end
