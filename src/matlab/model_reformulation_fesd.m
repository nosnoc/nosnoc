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
    settings.couple_across_stages = 0;
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
%% Stewart's representation of the sets R_i and discirimant functions g_i
g_ind_all = [ ];
c_all = [];
m_vec = [];

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
    % Check if all data is avilable and if dimensions match.
    if ~iscell(S)
        S = {S};
    end
    if length(S) ~= n_simplex
        error('Number of matrices S does not match number of subsystems (taken to be number of matrices F_i which collect the modes of every subsystem).')
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
    % Create Stewart's indicator functions g_ind_ii
    for ii = 1:n_simplex
        if size(S{ii},2) ~= length(c{ii})
            error('The matrix S and vector c do not have compatible dimension.');
        end
        % discrimnant functions
        eval(['g_ind_' num2str(ii) '= -S{ii}*c{ii};']);
        eval(['g_ind_all = [g_ind_all;' 'g_ind_' num2str(ii) '];']);
        eval(['c_all = [c_all; c{ii}];']);
    end
end

% index sets and dimensions for ubsystems
m_ind_vec = [cumsum(m_vec)-m_1+1]; % index ranges of the corresponding thetas and lambdas
m = sum(m_vec);
%% Parameters
if casadi_symbolic_mode == 'MX'
    sigma = MX.sym('sigma'); % homotopy parameter
else
    sigma = SX.sym('sigma');
end

% eval(['sigma = ' casadi_symbolic_mode '.sym(''sigma'');'])
n_param = 1;  % number of parameters,  we model it as control variables and merge them with simple equality constraints
p = [sigma];
n_p = 1;

%% Declare model variables and equations
% Algebraic variables defintion
n_theta = sum(m_vec); % number of modes
n_lambda = n_theta;
n_z = n_theta+n_lambda+n_simplex; % n_theta + n_lambda + n_mu

theta = [];
mu = [];
lambda = [];

E = []; % diagonal matrix with one vectors

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

%% Define algerbraic variables which arise from Stewart's reformulation of a PSS into a DCS
% symbolic variables
z = [theta;lambda;mu];
% inital guess for z
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

model.z0 = z0;
model.lbz = lbz;
model.ubz = ubz;


%% Reformulate the Filippov ODE into a DCS
f_x = zeros(n_x,1);
% rhs of ODE;
for ii = 1:n_simplex
    ii_str = num2str(ii);
    eval(['f_x = f_x + F_' ii_str '*theta_'  ii_str  ';']);

end

% basic algebraic equations and complementarty condtions of the DCS
% (Note that the cross complementarities are later defined when the discrete
% time variables for every IRK stage in the create_nlp_fesd function are defined.)

% g_ind_i - lambda_i + mu_i e_i = 0; for all i = 1,..., n_simplex
% lambda_i'*theta_i = 0; for all i = 1,..., n_simplex
% lambda_i >= 0;    for all i = 1,..., n_simplex
% theta_i >= 0;     for all i = 1,..., n_simplex

f_z = []; % collects standard algebraic equations 0 = g_i(x) - \lambda_i - e \mu_i
f_z_convex = []; % equation for the convex multiplers 1 = e' \theta
f_comp_residual = 0; % the orthogonality conditions diag(\theta) \lambda = 0.
for i = 1:n_simplex
    i_str = num2str(i);
    % Gradient of Lagrange Function of indicator LP
    eval(['f_z = [f_z; g_ind_' i_str '-lambda_'  i_str '+mu_'  i_str '*e_' i_str '];']);
    % Sum of all thethas equal to 1
    eval(['f_z_convex = [f_z_convex; ; e_' i_str '''*theta_'  i_str '-1];']);
    eval(['f_comp_residual = f_comp_residual + lambda_' i_str '''*theta_'  i_str ';']);
    %     end
end
% f_z = [f_z;f_z_convex];
g_lp = [f_z;f_z_convex];

%% MPCC Specific Considerations
% sum over all complementariteies
J_cc = f_comp_residual;  % (used in l1 penalties and for evaluation of resiudal)
% Point-wise
n_algebraic_constraints = n_theta+n_simplex;  % dim(g  - lambda - mu *e ) + dim( E theta) ;

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
    %     f_z_fun = Function('f_z_fun',{x,z,u},{f_z}); % old name
    g_lp_fun = Function('g_lp_fun',{x,z,u},{g_lp}); % lp kkt conditions without bilinear complementarity term (it is treated with the other c.c. conditions)
    f_J_cc = Function('f_J_cc',{x,z,u},{J_cc});
else
    f_x_fun = Function('f_x_fun',{x,z},{f_x,f_q});
    %     f_z_fun = Function('f_z_fun',{x,z},{f_z});
    g_lp_fun = Function('g_lp_fun',{x,z},{g_lp}); % lp kkt conditions without bilinear complementarity term (it is treated with the other c.c. conditions)
    f_J_cc = Function('f_J_cc',{x,z},{J_cc});
end

f_q_T_fun = Function('f_q_T',{x},{f_q_T});


%% Function for continious time output of algebraic variables.
% %% TODO: remove this after new treatmanet of boudnary is introduced;
% try
% [B,C,D,tau_root] = collocation_times_fesd(n_s,irk_scheme);
% tau_col = tau_root(2:end);
% eval(['tau = ' casadi_symbolic_mode '.sym(''tau'');'])
% eval(['h_scale = ' casadi_symbolic_mode '.sym(''h_scale'');']) % Rescale the intevral from [0 1] to [0 h_scale]; h_scale will usually be the step size
% eval(['y = ' casadi_symbolic_mode '.sym(''y'',n_s);'])
% eval(['y_vec = ' casadi_symbolic_mode '.sym(''y_vec'',n_s,n_theta);'])
% 
% lagrange_poly = 0;
% m_lagrange = [];
% for ii = 1:n_s
%     term_ii = y(ii);
%     m_temp = 1;
%     tau_scaled= h_scale*tau;
%     tau_col_scaled = h_scale*tau_col;
%     for jj = 1:n_s
%         if jj~=ii
%             term_ii = term_ii.*((tau_scaled-tau_col_scaled(jj))./(tau_col_scaled(ii)-tau_col_scaled(jj)));
%             m_temp  = m_temp*((1-tau_col(jj))./(tau_col(ii)-tau_col(jj)));
%         end
%     end
%     m_lagrange = [m_lagrange,m_temp];
%     lagrange_poly = lagrange_poly + term_ii;
% end
% lagrange_poly_fun = Function('lagrange_poly',{y,tau,h_scale},{lagrange_poly}); % note that this is a scalar function
% catch
% 
% end
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
model.f_z = f_z;
model.f_q_T = f_q_T;

model.f_x_fun = f_x_fun;
% model.f_z_fun = f_z_fun;
model.g_lp_fun = g_lp_fun;
model.f_q_T_fun = f_q_T_fun;

% model.f_z_cc = f_z_cc;
model.f_J_cc = f_J_cc;
model.g_ind_all_fun = g_ind_all_fun;
model.c_fun = c_fun;

% try
% model.lagrange_poly_fun = lagrange_poly_fun;
% model.m_lagrange = m_lagrange;
% catch
%     
% end
% Some Dimensions;
model.n_x = n_x;
model.n_z = n_z;
model.n_u = n_u;
model.n_p = n_p;
model.n_simplex = n_simplex;

model.z = z;
model.E = E;

model.m_vec = m_vec;
model.m_ind_vec = m_ind_vec;
model.n_theta = n_theta;
model.n_lambda = n_lambda;
model.n_algebraic_constraints = n_algebraic_constraints;


%% collect all dimensions in one sperate struct as it is needed by several other functions later.
dimensions.N_stages = N_stages;
dimensions.N_finite_elements = N_finite_elements;
dimensions.n_x = n_x;
dimensions.n_u = n_u;
dimensions.n_z = n_z;
dimensions.n_s = n_s;
dimensions.n_theta = n_theta;
dimensions.n_simplex = n_simplex;
dimensions.m_vec = m_vec;
dimensions.m_ind_vec = m_ind_vec;
model.dimensions = dimensions;
end
