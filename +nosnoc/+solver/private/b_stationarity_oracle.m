function solution = b_stationarity_oracle(nlp,problem_data,settings)
%% Implemented by Armin Nurkanovic, Anton Pozharskiy 2023
import casadi.*
%% extract inital data
%% check nlp input
if isfield(nlp,'x')
    x = vertcat(nlp.x{:});
    n_primal = length(x);
    % inital guess
    if isfield(problem_data,'x0')
        x_k = problem_data.x0;
    else
        x_k = ones(n_primal,1);
    end
else
    error('Please provide the vector of degrees of freedom x.')
end

%% check nlp input
if isfield(nlp,'f')
    f = nlp.f;
else
    f = 0;
    warning('No objective was provided.');
end

if isfield(nlp,'g')
    g = nlp.g;
    n_constraints = length(g);
else
    n_constraints = 0;
end

if isfield(nlp,'comp1') && isfield(nlp,'comp2')
    x1 = nlp.comp1;
    x2 = nlp.comp2;
    n_comp = size(x1,1);
    if length(x1)~= length(x2)
        error('The vector comp1 and comp2 must have the same size.')
    end
    % find index set
    if n_comp > 0
        ind_x1 = [];
        ind_x2 = [];
        tic
        ind_x1_fun = Function('ind_1',{x},{x.jacobian(x1)});
        [ind_x1,~] = find(sparse(ind_x1_fun(x_k)==1));
        ind_x2_fun = Function('ind_2',{x},{x.jacobian(x2)});
        [ind_x2,~] = find(sparse(ind_x2_fun(x_k))==1);

        settings.nlp_is_mpcc = 1;
    else
        settings.nlp_is_mpcc = 0;
        ind_x1 = [];
        ind_x2 = [];
    end
else
    settings.nlp_is_mpcc= 0;
end

n_primal_x0 = n_primal - 2*n_comp; % primal variables excluding the complementarity variables;
ind_x0 = 1:n_primal;
ind_x0 = ind_x0';
if settings.nlp_is_mpcc
    ind_x0([ind_x1,ind_x2]) = [];
end

x0 = x(ind_x0);



%% Parameters
if isfield(nlp,'p')
    p = nlp.p;
    if isfield(problem_data,'p')
        p_val = problem_data.p;
        if size(p_val,1)~=size(p,1)
            error('Length of p and its value must be the same.')
        end
    end
else
    p = [];
    p_val = [];
end

%% check problem data input and read
if any(problem_data.lbg > problem_data.ubg)
    error('For some component it holds that lbg > ubg.');
else
    lbg = problem_data.lbg;
    ubg = problem_data.ubg;
end
if any(problem_data.lbx > problem_data.ubx)
    error('For some component it holds that lbx > ubx.');
else
    lbx = problem_data.lbx;
    ubx = problem_data.ubx;
end

if length(lbx)~= length(x) || length(ubx) ~= length(x)
    error('lbx, ubx and x must be vectors of the same length.')
end

if length(lbg)~= length(g) || length(ubg) ~= length(g)
    error('lbg, ubg and g must be vectors of the same length.')
end

% update complementarity lowerbounds if the do not exist
if settings.nlp_is_mpcc
    lbx(ind_x1 == -inf) = 0;
    lbx(ind_x2 == -inf) = 0;
end


%% Split into equalites and inequalities
ind_g_eq = find(lbg == ubg);
ind_g_ineq = find(lbg ~= ubg);

ind_g_ineq_lb = find(lbg>-inf & lbg ~= ubg);
ind_g_ineq_ub = find(ubg<inf & lbg ~= ubg);

ind_x_lb = find(lbx>-inf);
ind_x_ub =  find(ubx<inf);

n_eq = length(ind_g_eq);
n_g_ineq_ub = length(ind_g_ineq_ub);
n_g_ineq_lb = length(ind_g_ineq_lb);

n_ubx = length(ind_x_ub);
n_lbx = length(ind_x_lb);

lbx_reduced = lbx(ind_x_lb);
ubx_reduced = ubx(ind_x_ub);

g_sym = g;
%% Primal variables
if settings.nlp_is_mpcc
    M = settings.BigM;
    Delta_TR_k = settings.Delta_TR_init;
    I_n_comp = eye(n_comp);
    B_comp_k = [zeros(n_comp,n_primal-2*n_comp), -eye(n_comp), zeros(n_comp), M*eye(n_comp);...
        zeros(n_comp,n_primal-2*n_comp), zeros(n_comp), -eye(n_comp), -M*eye(n_comp)];     % matrix for expressing complementarity constraints with integers
    if settings.tighten_bounds_in_lpcc
        B_comp_k = [zeros(2*n_comp,n_primal),[diag(min(ubx(ind_x1),M)); diag(-min(ubx(ind_x2),M))]];
    else
        B_comp_k = [zeros(2*n_comp,n_primal),[M*I_n_comp; -M*I_n_comp]];
    end

    for ii = 1:n_comp
        B_comp_k(ii,ind_x1(ii)) = -1;
        B_comp_k(n_comp+ii,ind_x2(ii)) = -1;
    end
    b_comp_k = zeros(2*n_comp,1); % to be updated at every iter

else
    B_comp_k  = zeros(0,n_primal);
    b_comp_k  = zeros(n_primal,0);
end

%% Constraint and objective jacobian
nabla_f = f.jacobian(x)';
% Zero order
g_eq = g(ind_g_eq)-lbg(ind_g_eq);                                          % g_eq = g - g_lb = 0
g_ineq_ub = ubg(ind_g_ineq_ub)-g(ind_g_ineq_ub);                           % g_ineq_ub = g_ub - g >= 0
g_ineq_lb = g(ind_g_ineq_lb)-lbg(ind_g_ineq_lb);                           % g_ineq_lb = g - g_lb >= 0

g_ineq = [ubg(ind_g_ineq_ub)-g(ind_g_ineq_ub);g(ind_g_ineq_lb)-lbg(ind_g_ineq_lb)];
n_ineq = size(g_ineq,1);

% first-order Constraint Jacobians
if n_eq > 0
    nabla_g = g_sym.jacobian(x);
    nabla_g_eq = g_eq.jacobian(x);
else
    nabla_g = [];
    nabla_g_eq = [];
end
if n_ineq > 0
    nabla_g_ineq_ub = g_ineq_ub.jacobian(x);
    nabla_g_ineq_lb = g_ineq_lb.jacobian(x);
else
    nabla_g_ineq_ub = [];
    nabla_g_ineq_lb = [];
end
nabla_g_ineq = [nabla_g_ineq_ub;nabla_g_ineq_lb];

%% CasADi functions of function evaluations and derivaties
% Zero order (objective and constraint function evaluations)
f_fun =  Function('f_fun',{x,p},{f});
g_eq_fun =  Function('g_eq_fun',{x,p},{g_eq});
g_ineq_ub_fun = Function('g_ineq_ub_fun',{x,p},{g_ineq_ub});
g_ineq_lb_fun = Function('g_ineq_lb_fun',{x,p},{g_ineq_lb});
g_ineq_fun = Function('g_ineq_lb_fun',{x,p},{g_ineq});

% First order (Gradients and Jacobian)
nabla_f_fun = Function('nabla_f_fun',{x,p},{nabla_f});
nabla_g_fun = Function('nabla_g_fun',{x,p},{nabla_g});
nabla_g_eq_fun = Function('nabla_g_eq_fun',{x,p},{nabla_g_eq});
nabla_g_ineq_ub_fun  = Function('nabla_g_ineq_ub_fun',{x,p},{nabla_g_ineq_ub});
nabla_g_ineq_lb_fun  = Function('nabla_g_ineq_lb_fun',{x,p},{nabla_g_ineq_lb});
nabla_g_ineq_fun  = Function('nabla_g_ineq_lb_fun',{x,p},{nabla_g_ineq});


%% Infeasiblity meausre in inf norm
% all standard constraints
h_eq = max(abs(g_eq));
h_ineq_ub = max(min(g_ineq_ub,0));
h_ineq_lb = max(min(g_ineq_lb,0));
h_ubx = max(min(ubx-x,0));
h_lbx = max(min(x-lbx,0));
% Summary
h_standard_constraints = max([h_eq;h_ineq_ub;h_ineq_lb;h_ubx;h_lbx]);
if n_comp > 0
    h_comp_constraints = max(abs(min(x1,x2)));
else
    h_comp_constraints  = 0;
end

%% CasADi Functions for constraints infeasiblity
h_standard_constraints_fun  = Function('h_standard_constraints_fun',{x,p},{h_standard_constraints});
h_comp_constraints_fun  = Function('h_comp_constraints_fun',{x,p},{h_comp_constraints});
h_total_fun = Function('h_comp_constraints_fun',{x,p},{max(h_comp_constraints,h_standard_constraints)});

%% Ininitalization
% Intalization step
d_k = zeros(n_primal,1);
d_lpcc = d_k;

if isfield(problem_data, 'y_lpcc')
    y_lpcc = problem_data.y_lpcc;
else
    y_lpcc = [];
end

%% Store some dimensions
dims.ind_x0  = ind_x0;
dims.ind_x1  = ind_x1;
dims.ind_x2  = ind_x2;
dims.n_primal = n_primal;
dims.n_comp = n_comp;
dims.n_eq = n_eq;
dims.n_ineq = n_ineq;

%% Stats
stopping_criterion_fullfiled = false;

k = 0;          % iteration counter

X = x_k; % save all iterates;
hessian_regularization_iter = [];
Delta_TR_iter = [];
f_iter = []; % objective
d_norm_iter = [];   % norm of step size
h_standard_iter = []; % standard infeasilibty
h_comp_iter = []; % comp infeasilibty

Delta_TR_k = settings.Delta_TR_init;
k_res  = 0; % number of trails of restoration phase

cpu_time_lpcc_iter = [];
cpu_time_qp_iter = [];
d_lpcc_list = [];
nabla_f_list = [];
%% Verbose

filter_message = ' ';
if settings.verbose_solver
    f_k = full(f_fun(x_k,p_val));
    h_standard_constraints_k = full(h_standard_constraints_fun(x_k,p_val));
    h_comp_constraints_k = full(h_comp_constraints_fun(x_k,p_val));
    fprintf('f_0 = %2.2e \n',f_k)
    fprintf('h_std_0 = %2.2e \n',h_standard_constraints_k )
    fprintf('h_comp_0 = %2.2e \n',h_comp_constraints_k )

    fprintf(' ------------------------------------------------------------------------\n');
    fprintf('D_tr\t\tf_lpcc\t\t||d||\t\tmax(abs(d))\n');
end
oracle_message  = 'Oracle failed to verify d=0 as locally optimal. \n';
oracle_status = 1;

%% Evaluate LPCC data
% Zero order
a_k  = full(g_eq_fun(x_k,p_val));
b_k  = full(g_ineq_fun(x_k,p_val));
f_k  = full(f_fun(x_k,p_val));

% First order
nabla_f_k  = full(nabla_f_fun(x_k,p_val));
nabla_f_list = [nabla_f_list,nabla_f_k];
A_k  = full(nabla_g_eq_fun(x_k,p_val));
B_k  = full(nabla_g_ineq_fun(x_k,p_val));

%% Main steps

Delta_TR_k_l = max(Delta_TR_k,settings.Delta_TR_init); % initalize TR
Delta_TR_iter  =  [Delta_TR_iter, Delta_TR_k];
l_k = 0; % iter counter
while l_k < settings.max_iter
    if settings.nlp_is_mpcc
        if settings.tighten_bounds_in_lpcc
            b_comp_k = [x_k(dims.ind_x1); x_k(dims.ind_x2)-min(ubx(dims.ind_x2),M)];
        else
            b_comp_k = [x_k(dims.ind_x1); x_k(dims.ind_x2)-M];
        end
    end
    
    Delta_TR_lpcc = Delta_TR_k_l;
    
    [f_lin_opt, d_lpcc_k, y_lpcc_k, lpcc_solution_exists, lpcc_cpu_time] = solve_lpcc(x_k,dims,...
        d_lpcc, y_lpcc, Delta_TR_lpcc,...
        lbx, ubx,...
        nabla_f_k, A_k, a_k, B_k, b_k, B_comp_k, b_comp_k);         % Solve LPCC
    cpu_time_lpcc_iter = [cpu_time_lpcc_iter, lpcc_cpu_time];

    if ~isempty(y_lpcc)
        activ_set_changes = sum(abs(y_lpcc_k-y_lpcc));
        %fprintf('activ set changes %d \n',activ_set_changes);
    end
    if settings.verbose_solver
        fprintf("%e\t%e\t%e\t%e\n", Delta_TR_lpcc, f_lin_opt, norm(d_lpcc_k), max(abs(d_lpcc_k)))
    end
    d_step_k = d_lpcc_k;
    d_lpcc_list = [d_lpcc_list, d_lpcc_k];
    l_k =  l_k + 1;
    if lpcc_solution_exists
        % Stopping criterion
        h_total_k = full(h_total_fun(x_k,p_val));
        % Note: we assume that the point is sufficiently feasible 
        if (h_total_k < settings.tol) && (abs(f_lin_opt) <= settings.tol || norm(nabla_f_k) < settings.tol)
            d_norm_iter = [d_norm_iter, norm(d_lpcc_k)];
            d_lpcc_k = d_lpcc_k*0;
        end

        if norm(d_lpcc_k) <= settings.tol && (h_total_k < settings.tol)
            d_norm_iter = [d_norm_iter, norm(d_lpcc_k)];
            oracle_message = 'B-stationary point found sucessfully! \n';
            oracle_status = 0;
            break;
        end
    end
    Delta_TR_k_l = Delta_TR_k_l/2;
end

% Final message
if settings.verbose_solver
    fprintf(oracle_message);
end

%% Extract Results
solution.x = x_k;
solution.f = f_k;
solution.d_lpcc = d_lpcc_list;
solution.nabla_f = nabla_f_list;

solution.oracle_status = oracle_status;
end
