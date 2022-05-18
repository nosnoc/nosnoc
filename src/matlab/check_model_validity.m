function [model] = check_model_validity(model,settings)
import casadi.*
unfold_struct(model,'caller');
%% Step-size
h = T/N_stages;
% nominal lengths of the finite elements for different control intevrals, every control interval might have a different number of finite elements.
N_finite_elements = N_finite_elements(:); % make a column vector of the input
if length(N_finite_elements) > N_stages
    N_finite_elements = N_finite_elements(1:N_stages);
    if print_level >=1
        fprintf('Info: Provided N_finite_elements had more antries then N_stages, the surplus of entries was removed. \n')
    end
end
if length(N_finite_elements) == 1
    N_finite_elements = N_finite_elements*ones(N_stages,1);
elseif length(N_finite_elements) > 1 && length(N_finite_elements) < N_stages
    N_finite_elements = N_finite_elements(:); % make sure it is a column vector
    N_finite_elements = [N_finite_elements;N_finite_elements(end)*ones(N_stages-length(N_finite_elements),1)];
end
h_k = h./N_finite_elements;

model.h = h;
model.h_k = h_k;
model.N_finite_elements = N_finite_elements;

%% Check is x

if exist('x')
    n_x = length(x);
    % check  lbx
    if exist('lbx')
        if length(lbx) ~= n_x
            error('The vector lbx, for the lower bounds of x has the wrong size.')
        end
    else
        lbx = -inf*ones(n_x,1);
    end
    % check ubx
    if exist('ubx')
        if length(ubx) ~= n_x
            error('The vector ubx, for the upper bounds of x has the wrong size.')
        end
    else
        ubx = inf*ones(n_x,1);
    end
else
    error('Please provide the state vector x, a CasADi symbolic variable.')
end
%% Check is u provided
if exist('u')
    n_u = length(u);
    % check  lbu
    if exist('lbu')
        if length(lbu) ~= n_u
            error('The vector lbu, for the lower bounds of u has the wrong size.')
        end
    else
        lbu = -inf*ones(n_u,1);
    end
    % check ubu
    if exist('ubu')
        if length(ubu) ~= n_u
            error('The vector ubu, for the upper bounds of u has the wrong size.')
        end
    else
        ubu = inf*ones(n_u,1);
    end
    % check u0
    if exist('u0')
        if length(u0) ~= n_u
            error('The vector u0, for the initial guess of u has the wrong size.')
        end
    else
        u0 = 0*ones(n_u,1);
    end
else
    n_u = 0;
    if print_level >=1
        fprintf('Info: No control vector u is provided. \n')
    end
    lbu = [];
    ubu = [];
end

%% Stage and terminal costs
if ~exist('f_q')
    if print_level >=1
        fprintf('Info: No stage cost is provided. \n')
    end
    eval(['f_q = ', casadi_symbolic_mode, '.zeros(1);'])
end
if exist('f_q_T')
    terminal_cost = 1;
else
    if print_level >=1
        fprintf('Info: No terminal cost is provided. \n')
    end
    eval(['f_q_T = ', casadi_symbolic_mode, '.zeros(1);'])
end

%% Inequality constraints
if exist('g_ineq')
    g_ineq_constraint  = 1;
    n_g_ineq = length(g_ineq);
    if exist('g_ineq_lb')
        if length(g_ineq_lb)~=n_g_ineq;
            error('The user provided vector g_ineq_lb has the wrong size.')
        end
    else
        g_ineq_lb = -inf*ones(n_g_ineq,1);
    end

    if exist('g_ineq_ub')
        if length(g_ineq_ub)~=n_g_ineq;
            error('The user provided vector g_ineq_ub has the wrong size.')
        end
    else
        g_ineq_ub =  0*ones(n_g_ineq,1);
    end
    g_ineq_fun  = Function('g_ineq_fun',{x,u},{g_ineq});
else
    n_g_ineq = 0;
    g_ineq_constraint  = 0;
    if print_level >=1
        fprintf('Info: No path constraints are provided. \n')
    end
end

%% Terminal constraints
if exist('g_terminal')
    terminal_constraint = 1;
    n_g_terminal = length(g_terminal);
    if exist('g_terminal_lb')
        if length(g_terminal_lb)~=n_g_terminal;
            error('The user provided vector g_terminal_lb has the wrong size.')
        end
    else
        g_terminal_lb = 0*ones(n_g_terminal,1);
    end

    if exist('g_terminal_ub')
        if length(g_terminal_ub)~=n_g_terminal;
            error('The user provided vector g_terminal_ub has the wrong size.')
        end
    else
        g_terminal_ub =  0*ones(n_g_terminal,1);
    end
    g_terminal_fun  = Function('g_terminal_fun',{x},{g_terminal});
else
    terminal_constraint = 0;
    n_g_terminal = 0;
    if print_level >=1
        fprintf('Info: No terminal constraints are provided. \n')
    end
end
%% Outputs
model.sigma = sigma;
model.p = p;
% 
model.lbx = lbx;
model.ubx = ubx;

model.lbu = lbu;
model.ubu = ubu;
if n_u > 0
    model.u0 = u0;
end

if g_ineq_constraint
    model.g_ineq_lb = g_ineq_lb;
    model.g_ineq_ub = g_ineq_ub;
    model.g_ineq_fun = g_ineq_fun;
    model.g_ineq_constraint = g_ineq_constraint;
end

if terminal_constraint
    model.g_terminal_lb = g_terminal_lb;
    model.g_terminal_ub = g_terminal_ub;
    model.terminal_constraint = terminal_constraint;
    model.g_terminal_fun = g_terminal_fun;
end

end