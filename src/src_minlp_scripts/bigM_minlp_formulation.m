%% big M MINLP formulation for bimodel system
import casadi.*

% Continuous time functions
f_A_fun = Function('f_A_fun', {x, u}, {f_A, L});
f_B_fun = Function('f_B_fun', {x, u}, {f_B});
f_q_T_fun = Function('f_q_T_fun', {x}, {f_q_T});
psi_fun = Function('psi_fun', {x}, {psi});
%% Start with an empty NLP
% binaires
n_y = 2;
M = 1e6;
eps = 1e-6;


y0 = zeros(n_y,1);
J = 0;

% "Lift" initial conditions
Xk = opti.variable(n_x);
opti.subject_to(Xk==x0);
opti.set_initial(Xk, x0);
ind_total = [1:n_x];
ind_x = [1:n_x];
ind_u = [];
ind_y = [];
ind_tf = [];

% Collect all states/controls
Xs = {Xk};
Us = {};
Ys = {};

discrete  = [];
discrete = [discrete; zeros(n_x,1)];
Y0 = [1;0];
if time_optimal_problem
    % New NLP variable for the control
    T_final = opti.variable();
    opti.subject_to(1e-2 <=T_final <=1e2);
    opti.set_initial(T_final, T_val);
    ind_tf = [ind_tf,ind_total(end)+1];
    ind_total = [ind_total,ind_total(end)+1];
    discrete = [discrete;0];
    T_val = 1;
    J = J + T_final;
    h = T/(N_stages*N_control_intevrlas*N_finite_elements);
else
    T_final = 1;
end

% Formulate the NLP
for k = 1:N_stages
    % define binary variables for stages / stage intervals
    if hybrid_dynamics
        Yk = opti.variable(n_y,1);
        ind_y = [ind_y,ind_total(end)+1:ind_total(end)+n_y];
        ind_total = [ind_total,ind_total(end)+1:ind_total(end)+n_y];

        YAk = Yk(1);
        YBk = Yk(2);
        Ys{end+1} = Yk;
        if k == 1
            opti.subject_to(Yk==Y0);
        end
        opti.set_initial(Yk, Y0);
        discrete = [discrete;ones(n_y,1)];
        % SOS1 constraint
        opti.subject_to(YAk + YBk == 1);
    end
    % define logical constraints regarding binary vairalbe implications
    for i = 1:N_control_intevrlas
        % New NLP variable for the control
        Uk = opti.variable(n_u);
        ind_u = [ind_u,ind_total(end)+1:ind_total(end)+n_u];
        ind_total = [ind_total,ind_total(end)+1:ind_total(end)+n_u];
        if g_u_constraint 
                       opti.subject_to(lb_g_u <= g_u_fun(Uk) <= ub_g_u);

        end

        Us{end+1} = Uk;
        opti.subject_to(lbu<=Uk<=ubu);
        opti.set_initial(Uk, u0);
        discrete = [discrete;0*ones(n_u,1)];
        for j = 1:N_finite_elements
            % Decision variables for helper states at each collocation point
            Xc = opti.variable(n_x, n_s);
            discrete = [discrete;zeros(n_x*n_s,1)];
            %             ind_x= [ind_x,ind_total(end)+1:ind_total(end)+n_x*n_s];
            ind_total = [ind_total,ind_total(end)+1:ind_total(end)+n_x*n_s];

            opti.subject_to(lbx <= Xc <=ubx);
            opti.set_initial(Xc, repmat(x0,1,n_s));
            % Evaluate ODE right-hand-side at all helper states
            [ode_A, quad] = f_A_fun(Xc, Uk);
            if hybrid_dynamics
                ode_B = f_B_fun(Xc, Uk);
                % Evaluate the transition functions
                [psi] = psi_fun(Xc);
            end

            % Add contribution to quadrature function
            J = J + quad*B*h;

            % Get interpolating points of collocation polynomial
            Z = [Xk Xc];

            % Get slope of interpolating polynomial (normalized)
            Pidot = Z*C;

            if hybrid_dynamics
                % Triggering of logical conditions
                opti.subject_to(-M*(YAk)<=psi<=(YBk)*M);
                % Match with ODE right-hand-side
                % Mode A
                opti.subject_to(-M*(1-YAk).*ones(n_s*n_x,1) <= Pidot(:) - T_final*h*ode_A(:) <= M*(1-YAk).*ones(n_s*n_x,1));
                % Mode B
                opti.subject_to(-M*(1-YBk).*ones(n_s*n_x,1) <= Pidot(:) - T_final*h*ode_B(:) <= M*(1-YBk).*ones(n_s*n_x,1));
            else
                % Match with ODE right-hand-side
                opti.subject_to(Pidot == T_final*h*ode_A);
            end

            % State at end of collocation interval
            Xk_end = Z*D;

            % New decision variable for state at end of interval
            Xk = opti.variable(n_x);
            ind_x= [ind_x,ind_total(end)+1:ind_total(end)+n_x];
            ind_total = [ind_total,ind_total(end)+1:ind_total(end)+n_x];

            discrete = [discrete;zeros(n_x,1)];
            Xs{end+1} = Xk;
            opti.subject_to(lbx <= Xk <= ubx);
            opti.set_initial(Xk, x0);
            % Continuity constraints
            opti.subject_to(Xk_end==Xk)
        end
    end
end
% terminal cost
J = J + f_q_T_fun(Xk);
n_terminal = length(x_goal);
if relax_terminal_constraint
    rho_terminal = 1e4;
    if 1
        s_ell1 = opti.variable(n_terminal);
        discrete = [discrete;zeros(n_terminal ,1)];
        opti.subject_to(-s_ell1  <= (Xk_end(1:n_terminal)-x_goal) <= s_ell1 );
        opti.set_initial(s_ell1 , 1e2);
        J = J+ rho_terminal*sum(s_ell1);
    else
        J = J+rho_terminal*(Xk_end(1:n_terminal)-x_goal)'*(Xk_end(1:n_terminal)-x_goal);
    end

else
    opti.subject_to(Xk_end(1:n_terminal) == x_goal);
end
opti.set_value(T, T_val)

opti.minimize(J);
%% collect variables
Xs = [Xs{:}];
Us = [Us{:}];
Ys = [Ys{:}];
