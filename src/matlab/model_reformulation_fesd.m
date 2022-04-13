%
%    This file is part of NOS-NOC.
%
%    NOS-NOC -- A software for NOnSmooth Numerical Optimal Control.
%    Copyright (C) 2022 Armin Nurkanovic, Moritz Diehl (ALU Freiburg).
%
%    NOS-NOC is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 3 of the License, or (at your option) any later version.
%
%    NOS-NOC is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with NOS-NOC; if not, write to the Free Software Foundation,
%    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%
function [model,settings] = model_reformulation_fesd(model,settings)
import casadi.*
%% Load settings and model details
unfold_struct(model,'caller');
unfold_struct(settings,'caller')

%% Sanity check of RK Schmes
if ~any(strcmp(list_of_all_rk_schemes,irk_scheme))
    if print_level >=1
        fprintf(['Info: The user provided RK scheme: ' irk_scheme ' is not supported, switching to Radau-IIA. \n'])
        fprintf(['See settings.list_of_all_rk_schemes in the default settings for an overview. \n'])
    end
    irk_scheme = 'Radau-IIA';
    settings.irk_scheme  = irk_scheme ;
end

if any(strcmp(list_of_all_rk_schemes(5:end),irk_scheme))
    if print_level >=1
        fprintf(['Info: The user provided RK scheme: ' irk_scheme ' is only avilabile in the differential representation.\n']);
    end
    irk_representation = 'differential';
    settings.irk_representation = irk_representation;
end

%% corrections for pss_mode
if isequal(settings.pss_mode,'stewart')  || isequal(settings.pss_mode,'stwrt') || isequal(settings.pss_mode,'indicator')
    settings.pss_mode = 'Stewart';
end

if isequal(settings.pss_mode,'step')  || isequal(settings.pss_mode,'stp') || isequal(settings.pss_mode,'Heaviside') ...
        || isequal(settings.pss_mode,'heaviside') || isequal(settings.pss_mode,'AP')
    settings.pss_mode = 'Step';
end

if ~isequal(settings.pss_mode,'Stewart') && ~isequal(settings.pss_mode,'Step')
    error('Please use for settings.pss_mode either ''Stewart''  or  ''Step''.' );
end

%% Some settings refinments
% update prin_level
if print_level <4
    settings.opts_ipopt.ipopt.print_level=0;
    settings.opts_ipopt.print_time=0;
    settings.opts_ipopt.ipopt.sb= 'yes';
elseif print_level == 4
    settings.opts_ipopt.ipopt.print_level=0;
    settings.opts_ipopt.print_time=1;
    settings.opts_ipopt.ipopt.sb= 'no';
else
    settings.opts_ipopt.ipopt.print_level = 5;
end

if settings.time_freezing
    settings.local_speed_of_time_variable = 1;
end

%% Determine is the SX or MX mode in CasADi used.
casadi_symbolic_mode = model.x(1).type_name();
settings.casadi_symbolic_mode  = casadi_symbolic_mode;
%% Step-size
h = T/N_stages;
% nominal lengths of the finite elements for different control intevrals,
% every control interval might have a different number of finite elements.
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

if n_u > 0
    settings.couple_across_stages = 1;
end

%% Stage and terminal costs
if ~exist('f_q')
    if print_level >=1
        fprintf('Info: No stage cost is provided. \n')
    end
    f_q = 0;
end
if exist('f_q_T')
    terminal_cost = 1;
else
    if print_level >=1
        fprintf('Info: No terminal cost is provided. \n')
    end
    f_q_T = 0;
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

%% Transforming a Piecewise smooth system into a DCS via Stewart's or the Step function approach
pss_mode = settings.pss_mode;


% Stewart's representation of the sets R_i and discirimant functions g_i
g_ind_all = [ ];
c_all = [];
m_vec = [];
n_c_vec = [];

if ~exist('F')
    error('Matrix F (or matrices F_i) with PSS modes not provided.');
else
    % check how many subsystems are present
    if iscell(F)
        n_simplex = length(F);
    else
        F = {F};
        n_simplex = 1;
    end
    % extract dimensions of subystems
    for ii = 1:n_simplex
        eval(['m_' num2str(ii) '= size(F{ii},2);']);
        eval(['m_vec = [m_vec m_' num2str(ii) '];'])
        eval(['F_' num2str(ii) '=F{ii};']);
    end
end

if ~exist('S')
    % if not the matrix S is provided, maybe the g_ind are avilable
    % directly?
    if isequal(pss_mode,'Stewart')
        if exist('g_ind')
            if ~iscell(g_ind)
                g_ind = {g_ind};
            end

            for ii = 1:n_simplex
                % discrimnant functions
                eval(['g_ind_' num2str(ii) '= g_ind{ii};']);
                eval(['g_ind_all = [g_ind_all;' 'g_ind_' num2str(ii) '];']);
                eval(['c_all = [c_all; 0];']);
            end
        else
            error(['Neither the sign matrix S nor the indicator functions g_ind for regions are provided. ' ...
                'Either provide the matrix S and the expression for c, or the expression for g_ind.']);
        end
    else
        error(['The user usses settings.pss_mode = ''Step'', but the sign matrix S is not provided. Please provide the matrix S and the expressions for c(x) (definfing the region boundaries).']);
    end
else
    % Check if all data is avilable and if dimensions match.
    if ~iscell(S)
        S = {S};
    end
    if length(S) ~= n_simplex
        error('Number of matrices S does not match number of subsystems. Note that the number of subsystems is taken to be number of matrices F_i which collect the modes of every subsystem.')
    end
    % Check constraint function c
    if ~exist('c')
        error('Expreesion for c, the constraint function for regions R_i is not provided.');
    else
        if ~iscell(c)
            c = {c};
        end
        if length(c) ~= n_simplex
            error('Number of different expressions for c does not match number of subsystems (taken to be number of matrices F_i which collect the modes of every subsystem).')
        end
    end

    % check are the matrices dense
    if isequal(pss_mode,'Stewart')
        for ii = 1:n_simplex
            if any(sum(abs(S{ii}),2)<size(S{ii},2))
                if n_simplex == 1
                    error('The matrix S is not dense. Either provide a dense matrix or use settings.mode = ''Step''.');
                else
                    error(['The matrix S{' num2str(ii) '} of the provided matrices is not dense. Either provide all dense matrices or use settings.mode = ''Step''.']);
                end
            end
        end
    end

    for ii = 1:n_simplex
        if size(S{ii},2) ~= length(c{ii})
            error('The matrix S and vector c do not have compatible dimension.');
        end

        % discrimnant functions
        switch pss_mode
            case 'Stewart'
                % Create Stewart's indicator functions g_ind_ii
                eval(['g_ind_' num2str(ii) '= -S{ii}*c{ii};']);
                eval(['g_ind_all = [g_ind_all;' 'g_ind_' num2str(ii) '];']);
            case 'Step'
                eval(['c_' num2str(ii) '= c{ii};']);
        end


        eval(['c_all = [c_all; c{ii}];']);
        % dimensions of c
        eval(['n_c_' num2str(ii) '= length(c{ii});']);
        n_c_vec  = [n_c_vec;length(c{ii})];
    end

end

% index sets and dimensions for ubsystems
m_ind_vec = [cumsum(m_vec)-m_1+1]; % index ranges of the corresponding thetas and lambdas
m = sum(m_vec);

%% Parameters
if isequal(casadi_symbolic_mode,'MX')
    sigma = MX.sym('sigma'); % homotopy parameter
else
    sigma = SX.sym('sigma');
end

% eval(['sigma = ' casadi_symbolic_mode '.sym(''sigma'');'])
n_param = 1;  % number of parameters,  we model it as control variables and merge them with simple equality constraints
p = [sigma];
n_p = 1;

%% Algebraic variables defintion
% dummy values for Stewart
theta = [];
mu = [];
lambda = [];
E = []; % diagonal matrix with one vectors

% dummy values for step
alpha  = [];
lambda_0 = [];
lambda_1 = [];
n_alpha = 0;
e_alpha = [];
n_beta = 0;
n_gamma = 0;
n_lambda_0 = 0;
n_lambda_1 = 0;

switch pss_mode
    case 'Stewart'
        % dimensions
        n_theta = sum(m_vec); % number of modes
        n_lambda = n_theta;
        n_f = n_theta;
        n_z = n_theta+n_lambda+n_simplex; % n_theta + n_lambda + n_mu


        % Define symbolic variables for algebraci equtions!
        for i = 1:n_simplex
            i_str = num2str(i);
            % define theta (convex multiplers of Filippov embedding)
            eval(['theta_' i_str '=' casadi_symbolic_mode '.sym(''theta_' i_str ''',m_' i_str ');'])
            eval(['theta = [theta; theta_' i_str '];'])
            % define mu_i (Lagrange multipler of e'theta =1;)
            eval(['mu_' i_str '= ' casadi_symbolic_mode '.sym(''mu_' i_str ''',1);'])
            eval(['mu = [mu; mu_' i_str '];'])
            % define lambda_i (Lagrange multipler of theta >= 0;)
            eval(['lambda_' i_str '= ' casadi_symbolic_mode '.sym(''lambda_' i_str ''',m_' i_str ');'])
            eval(['lambda = [lambda; lambda_' i_str '];'])
            % adefine ppropiate vector of ones
            eval(['e_' i_str '=ones(m_' i_str ' ,1);'])
            if n_simplex > 1
                eval(['E_row = zeros(m_' i_str ',n_simplex);'])
                eval(['E_row(:,' i_str ') = e_' i_str ';']);
                E = [E;E_row];
            end
        end
        if n_simplex == 1
            E = e_1;
        end
    case 'Step'
        n_alpha = sum(n_c_vec);
        n_f = sum(m_vec);
        n_lambda_0 = sum(n_c_vec);
        n_lambda_1 = sum(n_c_vec);
        % for creae_nlp_fesd
        n_theta = 2*n_alpha;
        n_lambda = n_lambda_0+n_lambda_1;
        % algebraic varaibles so far
        n_z = n_alpha+n_lambda_0+n_lambda_1;

        for i = 1:n_simplex
            i_str = num2str(i);
            % define alpha (selection of a set valued step function)
            eval(['alpha_' i_str '=' casadi_symbolic_mode '.sym(''alpha_' i_str ''',n_c_' i_str ');'])
            eval(['alpha = [alpha; alpha_' i_str '];'])
            % define lambda_0_i (Lagrange multipler of alpha >= 0;)
            eval(['lambda_0_' i_str '=' casadi_symbolic_mode '.sym(''lambda_0_' i_str ''',n_c_' i_str ');'])
            eval(['lambda_0= [lambda_0; lambda_0_' i_str '];'])
            % define lambda_1_i (Lagrange multipler of alpha <= 1;)
            eval(['lambda_1_' i_str '=' casadi_symbolic_mode '.sym(''lambda_1_' i_str ''',n_c_' i_str ');'])
            eval(['lambda_1= [lambda_1; lambda_1_' i_str '];'])
        end
        % adefine ppropiate vector of ones % for the kkt conditions of the LP
        e_alpha = ones(n_alpha,1);

        % Define already here lifting variables and functions
        beta = [];
        gamma = [];

        % Upsilo collects the vector for dotx = F(x)Upsilon, it is either multiaffine
        % terms or gamma from lifting
        pss_lift_step_functions = 0;
        if pss_lift_step_functions
            % do the lifting algortim for ever subystem
            % introduce varibles beta and gamma
            % introduce their dimensions
        else
            %             upsilon_i = {};
            % define the (multi)afine term for all subsystems
            for i = 1:n_simplex
                i_str = num2str(i);
                eval(['upsilon_' i_str ' = [];'])
                S_temp = S{i};
                for j = 1:size(S_temp,1)
                    upsilon_ij = 1;
                    for k = 1:size(S_temp,2)
                        % create multiafine term
                        if S_temp(j,k) ~=0
                            eval(['upsilon_ij = upsilon_ij * ( 0.5*(1-S_temp(j,k))+S_temp(j,k)*alpha_' i_str ' (k) ) ;'])
                        end
                    end
                    eval(['upsilon_' i_str ' = [upsilon_' i_str ';upsilon_ij ];'])
                end
                % go thorug rows of of every matrix and deffine the
                % appropate term
                % store them in the upsilon cell
            end
        end
        n_beta = length(beta);
        n_gamma = length(gamma);
        n_z = n_z + n_beta+n_gamma;
end


%% Define algerbraic variables which arise from Stewart's reformulation of a PSS into a DCS
switch pss_mode
    case 'Stewart'
        % symbolic variables
        z = [theta;lambda;mu];
        lbz = [0*ones(n_theta,1);0*ones(n_theta,1);-inf*ones(n_simplex,1)];
        ubz = [inf*ones(n_theta,1);inf*ones(n_theta,1);inf*ones(n_simplex,1)];
        % inital guess for z; % solve LP for guess;
        if lp_initalization
            [theta_guess,lambda_guess,mu_guess] = create_lp_based_guess(model);
        else
            theta_guess = initial_theta*ones(n_theta,1);
            lambda_guess = initial_lambda*ones(n_theta,1);
            mu_guess = initial_mu*ones(n_simplex,1);
        end
        z0 = [theta_guess;lambda_guess;mu_guess];
    case 'Step'
        z = [alpha;lambda_0;lambda_1;beta;gamma];
        lbz = [0*ones(n_alpha,1);0*ones(n_alpha,1);0*ones(n_alpha,1);-inf*ones(n_beta,1);-inf*ones(n_gamma,1)];
        ubz = [ones(n_alpha,1);inf*ones(n_alpha,1);inf*ones(n_alpha,1);inf*ones(n_beta,1);inf*ones(n_gamma,1)];

        alpha_guess = initial_alpha*ones(n_alpha,1);
        lambda_0_guess = initial_lambda_0*ones(n_alpha,1);
        lambda_1_guess = initial_lambda_1*ones(n_alpha,1);
        beta_guess = initial_beta*ones(n_beta,1);
        gamma_guess = initial_gamma*ones(n_gamma,1);
        % eval functios for gamma and beta?
        z0 = [alpha_guess;lambda_0_guess;lambda_1_guess;beta_guess;gamma_guess];
end

model.z0 = z0;
model.lbz = lbz;
model.ubz = ubz;

%% Reformulate the Filippov ODE into a DCS
f_x = zeros(n_x,1);
% rhs of ODE;

for ii = 1:n_simplex
    ii_str = num2str(ii);
    switch pss_mode
        case 'Stewart'
            eval(['f_x = f_x + F_' ii_str '*theta_'  ii_str  ';']);
        case 'Step'
            eval(['f_x = f_x + F_' ii_str '*upsilon_'  ii_str  ';']);
    end
end


g_z = []; % collects standard algebraic equations 0 = g_i(x) - \lambda_i - e \mu_i
g_z_convex = []; % equation for the convex multiplers 1 = e' \theta
g_lift = [];
f_comp_residual = 0; % the orthogonality conditions diag(\theta) \lambda = 0.



for i = 1:n_simplex
    i_str = num2str(i);
    switch pss_mode
        case 'Stewart'
            % basic algebraic equations and complementarty condtions of the DCS
            % (Note that the cross complementarities are later defined when the discrete
            % time variables for every IRK stage in the create_nlp_fesd function are defined.)
            % g_ind_i - lambda_i + mu_i e_i = 0; for all i = 1,..., n_simplex
            % lambda_i'*theta_i = 0; for all i = 1,..., n_simplex
            % lambda_i >= 0;    for all i = 1,..., n_simplex
            % theta_i >= 0;     for all i = 1,..., n_simplex
            % Gradient of Lagrange Function of indicator LP
            eval(['g_z = [g_z; g_ind_' i_str '-lambda_'  i_str '+mu_'  i_str '*e_' i_str '];']);
            eval(['g_z_convex = [g_z_convex; ; e_' i_str '''*theta_'  i_str '-1];']);
            eval(['f_comp_residual = f_comp_residual + lambda_' i_str '''*theta_'  i_str ';']);
        case 'Step'
            % c_i(x) - (lambda_1_i-lambda_0_i)  = 0; for all i = 1,..., n_simplex
            % lambda_0_i'*alpha_i  = 0; for all i = 1,..., n_simplex
            % lambda_1_i'*(e-alpha_i)  = 0; for all i = 1,..., n_simplex
            % lambda_0_i >= 0;    for all i = 1,..., n_simplex
            % lambda_1_i >= 0;    for all i = 1,..., n_simplex
            % alpha_i >= 0;     for all i = 1,..., n_simplex
            eval(['g_z = [g_z; c_' i_str '-lambda_1_'  i_str '+lambda_0_'  i_str '];']);
            % to do , greate g_lift
            eval(['f_comp_residual = f_comp_residual + lambda_0_' i_str '''*alpha_'  i_str '+ lambda_1_' i_str '''*(e_alpha-alpha_'  i_str ');']);
    end
end
g_lp = [g_z;g_z_convex;g_lift];
n_algebraic_constraints  = length(g_lp);
% n_algebraic_constraints = n_theta+n_simplex;  % dim(g  - lambda - mu *e ) + dim( E theta) ;

%% CasADi functions for indictaor and region constraint functions
% model equations
if n_u >0
    g_ind_all_fun = Function('g_ind_all_fun',{x,u},{g_ind_all});
    c_fun = Function('c_fun',{x,u},{c_all});
else
    g_ind_all_fun = Function('g_ind_all_fun',{x},{g_ind_all});
    c_fun = Function('c_fun',{x},{c_all});
end

if n_u >0
    f_x_fun = Function('f_x_fun',{x,z,u},{f_x,f_q});
    g_lp_fun = Function('g_lp_fun',{x,z,u},{g_lp}); % lp kkt conditions without bilinear complementarity term (it is treated with the other c.c. conditions)
else
    f_x_fun = Function('f_x_fun',{x,z},{f_x,f_q});
    g_lp_fun = Function('g_lp_fun',{x,z},{g_lp}); % lp kkt conditions without bilinear complementarity term (it is treated with the other c.c. conditions)
end

J_cc_fun = Function('J_cc_fun',{z},{f_comp_residual});
f_q_T_fun = Function('f_q_T',{x},{f_q_T});

%% Intigal guess for state derivatives at stage points
if isequal(irk_representation,'differential')
    if simple_v0_guess
        v0 = zeros(n_x,1);
    else
        if n_u>0
            [v0,~] = (f_x_fun(x0,z0,u0));
            v0 = full(v0);
        else
            [v0,~] = (f_x_fun(x0,z0));
            v0 = full(v0);
        end
    end
    model.v0 = v0;
end

%% Collect Outputs
model.sigma = sigma;
model.p = p;

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


model.f_x = f_x;
model.g_z = g_z;
model.g_lp = g_lp;
model.f_q_T = f_q_T;

model.f_x_fun = f_x_fun;
model.g_lp_fun = g_lp_fun;
model.f_q_T_fun = f_q_T_fun;

model.J_cc_fun = J_cc_fun;
model.g_ind_all_fun = g_ind_all_fun;
model.c_fun = c_fun;

% Model Dimensions;
model.n_x = n_x;
model.n_z = n_z;
model.n_u = n_u;
model.n_p = n_p;
model.n_simplex = n_simplex;

model.z = z;
model.E = E;
model.e_alpha = e_alpha;

model.m_vec = m_vec;
model.m_ind_vec = m_ind_vec;
model.n_theta = n_theta;
model.n_lambda = n_lambda;
model.n_algebraic_constraints = n_algebraic_constraints;

model.n_c_vec = n_c_vec;
model.n_alpha = n_alpha;
model.n_beta = n_beta;
model.n_gamma = n_gamma;
model.n_lambda_0 = n_lambda_0;
model.n_lambda_1 = n_lambda_1;


%% collect all dimensions in one sperate struct as it is needed by several other functions later.
dimensions.N_stages = N_stages;
dimensions.N_finite_elements = N_finite_elements;
dimensions.n_x = n_x;
dimensions.n_f = n_f;
dimensions.n_u = n_u;
dimensions.n_z = n_z;
dimensions.n_s = n_s;
dimensions.n_theta = n_theta;
dimensions.n_simplex = n_simplex;
dimensions.m_vec = m_vec;
dimensions.m_ind_vec = m_ind_vec;
dimensions.n_c_vec = n_c_vec;
dimensions.n_alpha = n_alpha;
dimensions.n_beta = n_beta;
dimensions.n_gamma = n_gamma;
dimensions.n_lambda_0 = n_lambda_0;
dimensions.n_lambda_1 = n_lambda_1;

model.dimensions = dimensions;
end
