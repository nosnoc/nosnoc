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
function [model,settings] = model_reformulation_nosnoc(model,settings)
import casadi.*
%% Load settings and model details
if ~settings.time_freezing_model_exists && settings.time_freezing && ~settings.time_freezing_hysteresis
    % check is the model generated if time freezing is used
    [model,settings] = time_freezing_reformulation(model,settings);
end

unfold_struct(model,'caller');
settings_bkp = settings;
unfold_struct(settings,'caller')
settings = settings_bkp;

%% Some settings refinments
% update prin_level
if print_level < 4
    settings.opts_ipopt.ipopt.print_level = 0;
    settings.opts_ipopt.print_time=0;
    settings.opts_ipopt.ipopt.sb= 'yes';
elseif print_level == 4
    settings.opts_ipopt.ipopt.print_level = 0;
    settings.opts_ipopt.print_time=1;
    settings.opts_ipopt.ipopt.sb= 'no';
else
    settings.opts_ipopt.ipopt.print_level = 5;
end

if settings.time_freezing
    settings.local_speed_of_time_variable = 1;
end

if isfield(model, 'alpha')
    settings.general_inclusion = 1;
else
    settings.general_inclusion = 0;
end

%% If different names are used...
if exist('N_stg','var')
    N_stages = N_stg;
    model.N_stages = N_stages;
end

if exist('N_FE','var')
    N_finite_elements = N_FE;
    model.N_finite_elements  = N_finite_elements;
end


%% Step-size
h = T/N_stages;
% nominal lengths of the finite elements for different control intevrals, every control interval might have a different number of finite elements.
N_finite_elements = N_finite_elements(:); % make a column vector of the input
if length(N_finite_elements) > N_stages
    N_finite_elements = N_finite_elements(1:N_stages);
    if print_level >=1
        fprintf('nosnoc: Provided N_finite_elements had more entries then N_stages, the surplus of entries was removed. \n')
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
%% Determine is the SX or MX mode in CasADi used.
casadi_symbolic_mode = class(model.x(1));
settings.casadi_symbolic_mode  = casadi_symbolic_mode;

%% Check is x provided
if isfield(model,'x')
    n_x = length(x);
    % check  lbx
    if isfield(model,'lbx')
        if length(model.lbx) ~= n_x
            error('nosnoc: The vector lbx, for the lower bounds of x has the wrong size.')
        end
    else
        lbx = -inf*ones(n_x,1);
    end
    % check ubx
    if isfield(model,'ubx')
        if length(model.ubx) ~= n_x
            error('nosnoc: The vector ubx, for the upper bounds of x has the wrong size.')
        end
    else
        ubx = inf*ones(n_x,1);
    end
else
    error('nosnoc: Please provide the state vector x, a CasADi symbolic variable.')
end
%% Check is u provided
if isfield(model,'u')
    n_u = length(model.u);
    % check  lbu
    if isfield(model,'lbu')
        if length(model.lbu) ~= n_u
            error('nosnoc: The vector lbu, for the lower bounds of u has the wrong size.')
        end
    else
        lbu = -inf*ones(n_u,1);
    end
    % check ubu
    if isfield(model,'ubu')
        if length(model.ubu) ~= n_u
            error('nosnoc: The vector ubu, for the upper bounds of u has the wrong size.')
        end
    else
        ubu = inf*ones(n_u,1);
    end
    % check u0
    if exist('u0')
        if length(u0) ~= n_u
            error('nosnoc: The vector u0, for the initial guess of u has the wrong size.')
        end
    else
        u0 = 0*ones(n_u,1);
    end
else
    u = define_casadi_symbolic(casadi_symbolic_mode,'',0);
    u0 = [];
    n_u = 0;
    if print_level >=1
        fprintf('nosnoc: No control vector u is provided. \n')
    end
    lbu = [];
    ubu = [];
end
%% Check if z is provided
if isfield(model,'z')
    n_z = length(z);

    if isfield(model,'z0')
        if length(model.z0) ~= n_z
            error('nosnoc: The vector z0, for the initial guess of z has the wrong size.')
        end
    else
        z0 = zeros(n_z, 1);
    end

    if isfield(model,'lbz')
        if length(model.lbz) ~= n_z
            error('nosnoc: The vector lbz, for the lower bound of z has the wrong size.')
        end
    else
        lbz = -inf*ones(n_z, 1);
    end

    if isfield(model,'ubz')
        if length(model.ubz) ~= n_z
            error('nosnoc: The vector ubz, for the lower bound of z has the wrong size.')
        end
    else
        ubz = inf*ones(n_z, 1);
    end
else
    n_z = 0;
    z0 = [];
    lbz = [];
    ubz = [];
    z = define_casadi_symbolic(casadi_symbolic_mode,'',0);
end
%% Global vars (i.e., variables that do not change with time)
if isfield(model,'v_global')
    n_v_global = length(model.v_global);
    if isfield(model,'v0_global')
        if length(model.v0_global) ~= n_v_global
            error('nosnoc: The vector v0_global, for the initial guess of v_global has the wrong size.')
        end
    else
        z0 = zeros(n_z, 1);
    end

    if isfield(model,'lbv_global')
        if length(model.lbv_global) ~= n_v_global
            error('nosnoc: The vector lbv_global, for the lower bound of v_global has the wrong size.')
        end
    else
        lbz = -inf*ones(n_z, 1);
    end

    if isfield(model,'ubv_global')
        if length(model.ubv_global) ~= n_v_global
            error('nosnoc: The vector ubv_global, for the upper bound of v_global has the wrong size.')
        end
    else
        ubz = inf*ones(n_z, 1);
    end
else
    n_v_global = 0;
    v_global = define_casadi_symbolic(casadi_symbolic_mode, '', 0);
    v0_global = [];
    lbv_global = [];
    ubv_global = [];
end

%% Parameters (time variable and that do not change with time)
if isfield(model,'p_global')
    n_p_global = size(model.p_global,1);
    if isfield(model,'p_global_val')
        if size(model.p_global_val,1) ~= n_p_global
            error('nosnoc: User provided p_global_val has the wrong size.')
        end
    else
        p_global_val = zeros(n_p_global,1);
    end
else
    n_p_global = 0;
    p_global = define_casadi_symbolic(casadi_symbolic_mode,'',0);
    p_global_val = [];
    if print_level >= 1
        fprintf('nosnoc: No global parameters given. \n')
    end
end

if isfield(model,'p_time_var')
    n_p_time_var = size(model.p_time_var, 1);
    if isfield(model,'p_time_var_val')
        if size(model.p_time_var_val) ~= [n_p_time_var, N_stages]
            error('nosnoc: User provided p_global_val has the wrong size.')
        end
    else
        p_time_var_val = zeros(n_p_time_var, N_stages);
    end

    p_time_var_stages = [];
    for ii=1:N_stages
        var_full = define_casadi_symbolic(casadi_symbolic_mode, ['p_time_var_' num2str(ii)], n_p_time_var);
        p_time_var_stages = horzcat(p_time_var_stages, var_full);
    end
else
    n_p_time_var = 0;
    p_time_var = define_casadi_symbolic(casadi_symbolic_mode,'',0);
    p_time_var_stages = define_casadi_symbolic(casadi_symbolic_mode,'', [0, N_stages]);
    p_time_var_val = double.empty(0,N_stages);
    if print_level >= 1
        fprintf('nosnoc: No time varying parameters given. \n')
    end
end

p = vertcat(p_global,p_time_var);

%% g_z: stage algebraic constraints
% TODO  long term: split up model_reformulation to allow f_alg to use the rest of stage Z
if isfield(model,'g_z')
    n_g_z = length(model.g_z);
else
    g_z = [];
    n_g_z = 0;
end
%% Stage and terminal costs check
if ~isfield(model,'f_q')
    if print_level >=1
        fprintf('nosnoc: No stage cost is provided. \n')
    end
    %     eval(['f_q = ', casadi_symbolic_mode, '.zeros(1);'])
    f_q = 0;
end

if isfield(model,'f_q_T')
    terminal_cost = 1;
else
    if print_level >=1
        fprintf('nosnoc: No terminal cost is provided. \n')
    end
    %     eval(['f_q_T = ', casadi_symbolic_mode, '.zeros(1);'])
    f_q_T = 0;
end
%% Least squares objective terms with variables references
if isfield(model,'lsq_x')
    if length(lsq_x)<3
        error('nosnoc: In lsq_x either the least squares function, the reference of the weight matrix are missing.')
    end
    if size(lsq_x{2},1)~=size(lsq_x{1})
        error('nosnoc: The dimensions of the least squares error term and weighting matrix for the differential states do not match.')
    end
    if size(lsq_x{1},1)~=size(lsq_x{3})
        error('nosnoc: The dimensions of the least squares error term and reference for the differential states do not match.')
    end

    n_x_ref_rows = size(lsq_x{2},1);
    n_x_ref_cols = size(lsq_x{2},2);
    if n_x_ref_cols == N_stages
        fprintf('nosnoc: the provided reference for the differential states is time variable. \n');
    elseif n_x_ref_cols == 1
        % replaciate
        fprintf('nosnoc: the provided reference for the differential states is constant over time. \n');
        lsq_x{2} = repmat(lsq_x{2},1,N_stages);
    else
        fprintf('nosnoc: The reference in lsq_x has to have a length of %d (if constant) or %d if time vriables. \n',1,N_stages)
        error('nosnoc: Please provide x_ref in lsq_x{1} with an appropaite size.')
    end
    x_ref_val = lsq_x{2};
    x_ref = define_casadi_symbolic(casadi_symbolic_mode,'x_ref',n_x_ref_rows);
    f_lsq_x = (lsq_x{1}-x_ref)'*lsq_x{3}*(lsq_x{1}-x_ref);
else
    x_ref = define_casadi_symbolic(casadi_symbolic_mode,'x_ref',1);
    f_lsq_x = 0;
    x_ref_val = zeros(1,N_stages);
end

% least square terms for control inputs
if isfield(model,'lsq_u')
    if length(model.lsq_u)<3
        error('nosnoc: In lsq_u either the least squares function, the reference of the weight matrix are missing.')
    end
    if size(model.lsq_u{2},1)~=size(model.lsq_u{1})
        error('nosnoc: The dimensions of the least squares error term and weighting matrix for the control input do not match.')
    end
    if size(model.lsq_u{1},1)~=size(model.lsq_u{3})
        error('nosnoc: The dimensions of the least squares error term and reference for the control input do not match.')
    end
    n_u_ref_rows = size(model.lsq_u{2},1);
    n_u_ref_cols = size(model.lsq_u{2},2);
    if n_u_ref_cols == N_stages
        fprintf('nosnoc: the provided reference for the control inputs is time variable. \n');
    elseif n_u_ref_cols == 1
        % replaciate
        fprintf('nosnoc: the provided reference for the control inputs is constant over time. \n');
        lsq_u{2} = repmat(lsq_u{2},1,N_stages);
    else
        fprintf('nosnoc: The reference in lsq_u has to have a length of %d (if constant) or %d if time vriables. \n',1,N_stages)
        error('nosnoc: Please provide u_ref in lsq_u{2} with an appropaite size.')
    end
    u_ref_val = lsq_u{2};
    u_ref = define_casadi_symbolic(casadi_symbolic_mode,'u_ref',n_u_ref_rows);
    f_lsq_u = (lsq_u{1}-u_ref)'*lsq_u{3}*(lsq_u{1}-u_ref);
else
    u_ref = define_casadi_symbolic(casadi_symbolic_mode,'u_ref',1);
    f_lsq_u = 0;
    u_ref_val = zeros(1,N_stages);
end


% least square terms for control inputs
if isfield(model,'lsq_T')
    % sanity chkecs on the input
    if length(model.lsq_T)<3
        error('nosnoc: In lsq_u either the least squares function, the reference of the weight matrix are missing.')
    end
    if size(model.lsq_T{2},1)~=size(model.lsq_T{1})
        error('nosnoc: The dimensions of the least squares error term and weighting matrix for the terminal cost do not match.')
    end
    if size(model.lsq_T{1},1)~=size(model.lsq_T{3})
        error('nosnoc: The dimensions of the least squares error term and reference for the terminal cost do not match.')
    end

    n_x_T_rows = size(lsq_T{2},1);
    n_x_T_cols = size(lsq_T{2},2);
    if n_x_T_cols == 1
        fprintf('nosnoc: the provided reference for the terminal cost is ok. \n');
    else
        fprintf('nosnoc: The reference in lsq_T has to be a vector of length %d. \n',length(lsq_T{1}));
        error('nosnoc: Please provide a reference vector in lsq_T{2} with an appropaite size.')
    end
    x_ref_end_val = lsq_T{2};
    x_ref_end = define_casadi_symbolic(casadi_symbolic_mode,'x_ref_end',n_x_T_rows);
    f_lsq_T = (lsq_T{1}-x_ref_end)'*lsq_T{3}*(lsq_T{1}-x_ref_end);
else
    x_ref_end  = define_casadi_symbolic(casadi_symbolic_mode,'x_ref_end',1);
    f_lsq_T = 0;
    x_ref_end_val = 0;
end

%% Inequality constraints check
if isfield(model,'g_path')
    g_path_constraint  = 1;
    n_g_path = length(model.g_path);
    if isfield(model,'g_path_lb')
        if length(model.g_path_lb)~=n_g_path
            error('nosnoc: The provided vector g_path_lb has the wrong size.');
        end
    else
        g_path_lb = -inf*ones(n_g_path,1);
    end

    %% Check path complementarity constraints
    g_comp_path_constraint  = 0;
    if isfield(model,'g_comp_path')
        g_comp_path_constraint  = 1;
        if size(g_comp_path, 2)
            error('g_comp_path must be of size (m, 2)')
        end
    else
        g_comp_path_constraint = 0;
        if print_level >=1
            fprintf('nosnoc: No path complementarity constraints are provided. \n')
        end
    else
        g_comp_path_ub =  0*ones(n_g_comp_path,1);
    end
    g_comp_path_fun  = Function('g_comp_path_fun',{x,u,p,v_global},{g_comp_path});
else
    n_g_comp_path = 0;
    g_comp_path_constraint  = 0;
    if print_level >=1
        fprintf('nosnoc: No path complementarity constraints are provided. \n')
    end
end
%% Terminal constraints
if isfield(model,'g_terminal')
    terminal_constraint = 1;
    n_g_terminal = length(model.g_terminal);
    if exist('g_terminal_lb')
        if length(g_terminal_lb)~=n_g_terminal
            error('nosnoc: The provided vector g_terminal_lb has the wrong size.')
        end
    else
        g_terminal_lb = 0*ones(n_g_terminal,1);
    end

    if isfield(model,'g_terminal_ub')
        if length(g_terminal_ub)~=n_g_terminal
            error('nosnoc: The provided vector g_terminal_ub has the wrong size.')
        end
    else
        g_terminal_ub =  0*ones(n_g_terminal,1);
    end
    g_terminal_fun  = Function('g_terminal_fun',{x,p_global,v_global},{g_terminal});
else
    terminal_constraint = 0;
    n_g_terminal = 0;
    if print_level >=1
        fprintf('nosnoc: No terminal constraints are provided. \n')
    end
end

%% Check are the inputs specific to a DCS mode correct and make some modifications
dcs_mode = settings.dcs_mode;
g_Stewart = {};
g_ind_vec = [];
c_all = [];
m_vec = [];
n_c_sys = [];

if isequal(dcs_mode,'CLS')
    % TODO: there is some repetition to the time_freezing check, this should be unified!!!!
    % Check existence of relevant functions
    n_sys = 1; % always one subystem in CLS (only loops over n_contacts later)
    if ~isfield(model,'f_c')
        error('nosnoc: Please provide the gap functions model.f_c.')
    end
    n_contacts = length(model.f_c);

    % coefficient of friction checks
    if isfield(model,'mu')
        if length(model.mu) ~= 1 || length(model.mu) ~= n_contacts
            error('The length of model.mu has to be one or match the length of model.f_c')
        end
        if length(model.mu) == 1
            model.mu = model.mu*ones(n_contacts,1);
        end

        if any(model.mu > 0)
            friction_exists = 1;
        else
            friction_exists = 0;
        end
    else
        model.mu = zeros(n_contacts,1);
        fprintf('nosnoc: Coefficients of friction mu not provided, setting it to zero for all contacts. \n')
    end
    if any(mu<0)
        error('nosnoc: The coefficients of friction mu should be nonnegative.')
    end

    % coefficent of restiution check
    if ~isfield(model,'e')
        error('nosnoc:  Please provide a coefficient of restitution via model.e')
    else
        if length(model.e) ~= 1 || length(model.e) ~= n_contacts
            error('The length of model.e has to be one or match the length of model.f_c')
        end
        if length(model.e) == 1
            model.e = model.e*ones(n_contacts,1);
        end
    end
    if any(abs(1-e)>1) || any(e<0)
        error('nosnoc: the coefficient of restitution e should be in [0,1].')
    end

    % dimensions and state space split
    casadi_symbolic_mode = model.x(1).type_name();
    if mod(size(x,1),2)
        n_x = size(x,1);
        n_q = (n_x-1)/2;
    else
        n_x = size(x,1);
        n_q = n_x/2;
    end
    if ~isfield(model,'q') && ~isfield(model,'v')
        q = x(1:n_q);
        v = x(n_q+1:2*n_q);
    end

    if ~isfield(model,'f_v')
        error('nosnoc: the function f_v (collecting all generalized forces), in M(q) = dv/dt =  f_v(q,v,u) + J_n\lambda_n +J_t\lambda_t ~ is not provided in model.');
    end

    % Check intertia matrix
    if ~isfield(model,'M')
        fprintf('nosnoc: Inertia matrix M is not provided, set to default value: M = eye(n_q). \n')
        M = eye(n_q);
        invM = inv(M);
    else
        invM = inv(M);
    end

    %  Normal Contact Jacobian
    if isfield(model,'J_normal')
        J_normal = model.J_normal;
        J_normal_exists = 1;
    else
        J_normal_exists = 0;
    end

    if J_normal_exists
        if size(J_normal,1)~=n_q && size(J_normal,2)~=n_contacts
            fprintf('nosnoc: J_normal should be %d x %d matrix.\n',n_q,n_contacts);
            error('nosnoc: J_normal has the wrong size.')
        end
        J_normal_exists = 1;
    else
        J_normal = f_c.jacobian(q)';
        fprintf('nosnoc: normal contact Jacobian not provided, but it is computed from the gap functions.\n');
        J_normal_exists = 1;
    end

    if is_zero(J_normal)
        error('nosnoc: The normal vector should have at least one non-zero entry.')
    end

    % Tangent Contact Jacobian
    if friction_exists
        if isequal(friction_model,'Conic')
            if isfield(model,'J_tangent')
                J_tangent = model.J_tangent;
                if size(J_tangent,1)~=n_q
                    error('nosnoc: J_tangent has the wrong size.')
                end
            else
                error('nosnoc: please provide the tangent Jacobian in model.J_tangent.')
            end
        end

        if isequal(friction_model,'Polyhedral')
            if ~isfield(model,'D_tangent')
                error('nosnoc: please provide the polyhedral tangent Jacobian in model.D_tangent, e.g., using the conic tangent Jacobian model.J_tangent: D_tangent = [J_tangent(q_0),-J_tangent(q_0)].')
            end
        end
    end
    % Dimension of tangens
    if friction_exists
        if isequal(friction_model,'Polyhedral')
            n_t =  length(D_tangent)/n_contacts; % number of tanget multipliers for a single contactl
        elseif isequal(friction_model,'Conic')
            n_t =  size(J_tangent,2)/n_contacts; % number of tanget multipliers for a single contactl
        else
            % TODO: @Anton, is this check needed at all, since we have the options class now.
            error('Pick friction_model Conic or Polyhedral')% position in symbolic vector z
        end
        n_tangents = n_t*n_contacts; % number tangent forces for all multpliers
    else
        n_tangents = 0;
    end
end
% TODO: the time freezing reformulation could be carried out at this point
% insetad at the begining and then go through the step/stewart checks.
if isequal(dcs_mode,'Step') || isequal(dcs_mode,'Stewart')
    if ~exist('F')
        % Don't need F
        if ~settings.general_inclusion
            error('nosnoc: Matrix F (or matrices F_i) with PSS modes not provided.');
        else
            % TODO Implement more subsystems.
            n_sys = 1;
            m_vec = [size(f_x,1)];
        end
    else
        % check how many subsystems are present
        if iscell(F)
            n_sys = length(F);
        else
            F = {F};
            n_sys = 1;
        end
        % extract dimensions of subystems
        for ii = 1:n_sys
            m_temp = size(F{ii},2);
            m_vec  = [m_vec m_temp];
        end
    end

    if ~exist('S')
        % if we are using general inclusions we dont need S.
        if ~settings.general_inclusion
            % if the matrix S is not provided, maybe the g_ind are available
            % directly?
            if isequal(dcs_mode,'Stewart')
                if exist('g_ind')
                    if ~iscell(g_ind)
                        g_ind = {g_ind};
                    end

                    for ii = 1:n_sys
                        % discriminant functions
                        g_ind_vec =  [g_ind_vec;g_ind{ii};];
                        g_Stewart{ii} = g_ind{ii};
                        c_all = [c_all; zeros(1,casadi_symbolic_mode)];
                    end
                else
                    error(['Neither the sign matrix S nor the indicator functions g_ind for regions are provided. ' ...
                        'Either provide the matrix S and the expression for c, or the expression for g_ind.']);
                end
            else
                error(['The user uses settings.dcs_mode = ''Step'', but the sign matrix S is not provided. Please provide the matrix S and the expressions for c(x) (definfing the region boundaries).']);
            end
        else
            if ~exist('c')
                error('nosnoc: Expresion for c, the constraint function for regions R_i is not provided.');
            else
                if ~iscell(c)
                    c = {c};
                end
                if length(c) ~= n_sys
                    error('nosnoc: Number of different expressions for c does not match number of subsystems.')
                end
                for ii = 1:n_sys
                    c_all = [c_all; c{ii}];
                    n_c{ii} = length(c{ii});
                    n_c_sys  = [n_c_sys;length(c{ii})];
                end

            end
        end
    else
        % Check if all data is avilable and if dimensions match.
        if ~iscell(S)
            S = {S};
        end
        if length(S) ~= n_sys
            error('nosnoc: Number of matrices S does not match number of subsystems. Note that the number of subsystems is taken to be number of matrices F_i which collect the modes of every subsystem.')
        end
        % Check constraint function c
        if ~exist('c')
            error('nosnoc: Expresion for c, the constraint function for regions R_i is not provided.');
        else
            if ~iscell(c)
                c = {c};
            end
            if length(c) ~= n_sys
                error('nosnoc: Number of different expressions for c does not match number of subsystems (taken to be number of matrices F_i which collect the modes of every subsystem).')
            end
        end

        % check are the matrices dense
        if isequal(dcs_mode,'Stewart')
            for ii = 1:n_sys
                if any(sum(abs(S{ii}),2)<size(S{ii},2))
                    if n_sys == 1
                        error('nosnoc: The matrix S is not dense. Either provide a dense matrix or use settings.mode = ''Step''.');
                    else
                        error(['The matrix S{' num2str(ii) '} of the provided matrices is not dense. Either provide all dense matrices or use settings.mode = ''Step''.']);
                    end
                end
            end
        end

        for ii = 1:n_sys
            if size(S{ii},2) ~= length(c{ii})
                error('nosnoc: The matrix S and vector c do not have compatible dimension.');
            end

            % discrimnant functions
            switch dcs_mode
                case 'Stewart'
                    % Create Stewart's indicator functions g_ind_ii
                    g_Stewart{ii} = -S{ii}*c{ii};
                    g_ind_vec = [g_ind_vec ;-S{ii}*c{ii}];
                case 'Step'
                    %eval(['c_' num2str(ii) '= c{ii};']);
            end
            % dimensions of c
            c_all = [c_all; c{ii}];
            n_c{ii} = length(c{ii});
            n_c_sys  = [n_c_sys;length(c{ii})];
        end

    end
    % index sets and dimensions for ubsystems
    m_ind_vec = 1;
    if isequal(dcs_mode,'Step')
        % double the size of the vectors, since alpha,1-alpha treated at same
        % time;
        m_vec = sum(n_c_sys)*2;
    end
    for ii = 1:length(m_vec)-1
        m_ind_vec  = [m_ind_vec,m_ind_vec(end)+m_vec(ii)];
    end
    % m_ind_vec = [cumsum(m_vec)-m_vec(1)+1]; % index ranges of the corresponding thetas and lambdas
    m = sum(m_vec);

    if isempty(n_c_sys)
        n_c_sys = 0;
    end

    if max(n_c_sys) < 2 && isequal(dcs_mode,'Step')
        pss_lift_step_functions = 0;
        if print_level >=1
            fprintf('nosnoc: settings.pss_lift_step_functions set to 0, as are step fucntion selections are already entering the ODE linearly.\n')
        end
    end

    if ~settings.general_inclusion
        n_f_sys = arrayfun(@(sys) size(F{sys},2),1:n_sys);
    else
        n_f_sys = [size(f_x,1)];
    end
end

%% TODO %%%%%%%%%%%%%%
%%%%%% HERE ENDS FIRST FUNCTION THAT JUST CHECKS THE CORECTNESS OF THE MODEL INPUT
%% TODO %%%%%%%%%%%%%%

%% Algebraic variables defintion
% Dummy variables for Stewart's representation
theta = [];
mu = [];
lambda = [];
% structs storing all vectors of every subsystem (they might have different dimensions)
theta_all = {};
lambda_all = {};
mu_all = {};
e_ones_all = {};

% dummy values for Step representation
if ~settings.general_inclusion
    alpha = [];
else
    alpha = vertcat(model.alpha{:});
end
lambda_n = [];
lambda_p = [];
alpha_all  = {};
lambda_0_all = {};
lambda_1_all = {};
theta_step_all = {};

% dummy values for CLS representation
lambda_normal = [];  % normal contact force
lambda_tangent = [];
gamma_d = [];
gamma = [];
beta = [];
p_vn = [];
n_vn = [];
p_vt = [];
n_vt = [];
alpha_vt = [];
n_gamma_d = 0;
n_gamma = 0;
% TODO: Add different flags to steer lifting
g_lift_gap = []; % lift gap functions f_c(q) = y;
g_lift_friction1 = []; % lift friction expression: Conic: beta = (mu \lambda_n)^2 - |\ \lambda_t \|_2^2 / Polyhedral: beta_d = mu \lambda_n - e^top \lambda_t
g_lift_friction2 = []; % lift stationarity condition in Polyhedral: detla_d = D_t(q)^\top v + \gamma_d e (becasue D(q) is nonlinear)

% dimensions
n_alpha = 0;
e_alpha = [];
n_beta = 0;
n_theta_step = 0;
n_lambda_n = 0;
n_lambda_p = 0;
g_lift_beta = [];
g_lift_theta_step  =[];
switch dcs_mode
    case 'Stewart'
        % dimensions
        n_theta = sum(m_vec); % number of modes
        n_lambda = n_theta;
        n_f = n_theta;
        n_z_all = n_theta+n_lambda+n_sys; % n_theta + n_lambda + n_mu
        % Define symbolic variables for algebraic equtions.
        for ii = 1:n_sys
            ii_str = num2str(ii);
            % define theta (Filippov multiplers)
            theta_temp = define_casadi_symbolic(casadi_symbolic_mode,['theta_' ii_str],m_vec(ii));
            theta = [theta;theta_temp];
            theta_all{ii} = theta_temp;
            % define mu_i (Lagrange multipler of e'theta =1;)
            mu_temp = define_casadi_symbolic(casadi_symbolic_mode,['mu_' ii_str],1);
            mu = [mu;mu_temp];
            mu_all{ii} = mu_temp;
            % define lambda_i (Lagrange multipler of theta >= 0;)
            lambda_temp = define_casadi_symbolic(casadi_symbolic_mode,['lambda_' ii_str],m_vec(ii));
            lambda = [lambda;lambda_temp];
            lambda_all{ii} = lambda_temp;
            % define appropiate vector of ones (the struct below stores them for every mode)
            e_ones_all{ii} = ones(m_vec(ii),1);
        end
    case 'Step'
        n_alpha = sum(n_c_sys);
        n_f = sum(m_vec);
        n_lambda_n = sum(n_c_sys);
        n_lambda_p = sum(n_c_sys);
        % for creae_nlp_fesd
        n_theta = 2*n_alpha;
        n_lambda = n_lambda_n+n_lambda_p;
        % algebraic varaibles so far
        n_z_all = n_alpha+n_lambda_n+n_lambda_p;
        for ii = 1:n_sys
            ii_str = num2str(ii);
            % define alpha (selection of a set valued step function)
            if ~settings.general_inclusion
                alpha_temp = define_casadi_symbolic(casadi_symbolic_mode,['alpha_' ii_str],n_c_sys(ii));
                alpha = [alpha;alpha_temp];
                alpha_all{ii} = alpha_temp;
            else
                % TODO this needs to change if subsystems.
                alpha_all{ii} = alpha{ii};
            end
            % define lambda_0_i (Lagrange multipler of alpha >= 0;)
            lambda_0_temp = define_casadi_symbolic(casadi_symbolic_mode,['lambda_0_' ii_str],n_c_sys(ii));
            lambda_n = [lambda_n;lambda_0_temp];
            lambda_0_all{ii} = lambda_0_temp;
            % define lambda_1_i (Lagrange multipler of alpha <= 1;)
            lambda_1_temp = define_casadi_symbolic(casadi_symbolic_mode,['lambda_1_' ii_str],n_c_sys(ii));
            lambda_p = [lambda_p;lambda_1_temp];
            lambda_1_all{ii} = lambda_1_temp;
        end
        % define appropiate vector of ones % for the kkt conditions of the LP
        e_alpha = ones(n_alpha,1);
        % Define already here lifting variables and functions
        % TODO allow for custom beta lifting
        beta = [];
        theta_step = [];
        % Theta collects the vector for dot_x = F(x)Theta ,
        % terms or theta_step from lifting;
        if ~settings.general_inclusion
            for ii = 1:n_sys
                theta_temp = [];
                ii_str = num2str(ii);
                S_temp = S{ii};
                if pss_lift_step_functions
                    % TODO implement automatic lifting
                else
                    if ~settings.time_freezing_inelastic
                        for j = 1:size(S_temp,1)
                            alpha_ij = 1;
                            for k = 1:size(S_temp,2)
                                % create multiafine term
                                if S_temp(j,k) ~=0
                                    alpha_ij = alpha_ij * (0.5*(1-S_temp(j,k))+S_temp(j,k)*alpha_all{ii}(k) ) ;
                                end
                            end
                            theta_temp = [theta_temp;alpha_ij];
                        end
                    end
                end
                theta_step_all{ii} = theta_temp;
            end
        end

        %% time-freezing inelastic impacts (exploit structure with taiolored formulae)
        if settings.time_freezing_inelastic
            % theta_step are the lifting variables that enter the ODE r.h.s.
            if ~settings.nonsmooth_switching_fun
                alpha_q = alpha(1:n_contacts);
                alpha_v_normal = alpha(n_contacts+1:2*n_contacts);
                if friction_exists
                    alpha_v_tangent = alpha(2*n_contacts+1:end);
                end
            else
                alpha_qv = alpha(1:n_contacts);
                if friction_exists
                    alpha_v_tangent = alpha(n_contacts+1:end);
                end
            end

            theta_step = define_casadi_symbolic(casadi_symbolic_mode,'theta_step',n_aux+1);
            theta_step_expr = eval([casadi_symbolic_mode '.zeros(' num2str(n_aux) '+1,1);']); % expressions for Filippov multipliers via alpha (and possibley beta).

            % lifting variables
            % empty expressions for initalization
            beta_bilinear_ode = []; % for lifting bilinear terms in free flight dynamics multiplire
            beta_bilinear_aux = []; % for lifiting bilinear terms appearing in aux. dynamics mutiplieres
            beta_prod = []; % for lifting the multi affine term defineng the overall free flight dynamics multpliers
            % expressions for lifting
            beta_bilinear_ode_expr = [];
            beta_bilinear_aux_expr = [];
            beta_prod_expr = [];
            beta_prod_expr_guess = []; % extra expresion to make depend only on alpha (the one above depens on both and alpha and beta) - needed for eval. of inital guess

            if pss_lift_step_functions
                % lift bilinear terms in product terms for free flight ode % (alpha_q*alpha_v)
                if ~nonsmooth_switching_fun
                    beta_bilinear_ode = define_casadi_symbolic(casadi_symbolic_mode,'beta_bilinear_ode',n_contacts);
                    beta_bilinear_ode_expr = eval([casadi_symbolic_mode '.zeros(' num2str(n_contacts) ',1);']);
                    if friction_exists
                        % lift bilinear terms defining aux dynamics (1-alpha_q)*(1-alpha_v)
                        beta_bilinear_aux = define_casadi_symbolic(casadi_symbolic_mode,'beta_bilinear_aux',n_contacts);
                        beta_bilinear_aux_expr = eval([casadi_symbolic_mode '.zeros(' num2str(n_contacts) ',1);']);
                    end
                end
                if n_contacts > 2
                    beta_prod = define_casadi_symbolic(casadi_symbolic_mode,'beta_bilinear',n_contacts-2);
                    beta_prod_expr = eval([casadi_symbolic_mode '.zeros(' num2str(n_contacts) '-2,1);']);
                    beta_prod_expr_guess = eval([casadi_symbolic_mode '.zeros(' num2str(n_contacts) '-2,1);']);
                end
            end
            beta = [beta_bilinear_ode;
                beta_bilinear_aux;
                beta_prod];
            % expresions for theta's and lifting
            %% Filippov multipliers
            alpha_ode = 1; % initalized product for free flight multiplier
            if ~pss_lift_step_functions
                for ii = 1:n_contacts
                    if nonsmooth_switching_fun
                        alpha_ode = alpha_ode*alpha_qv(ii);
                        if friction_exists
                            theta_step_expr(ii+1) = (1-alpha_qv(ii))*(alpha_v_tangent(ii));
                            theta_step_expr(1+n_contacts+ii) = (1-alpha_qv(ii))*(1-alpha_v_tangent(ii));
                        else
                            theta_step_expr(ii+1)=(1-alpha_qv(ii));
                        end
                    else
                        alpha_ode = alpha_ode*(alpha_q(ii)+alpha_v_normal(ii)-alpha_q(ii)*alpha_v_normal(ii));
                        if friction_exists
                            theta_step_expr(ii+1) = (1-alpha_q(ii))*(1-alpha_v_normal(ii))*(alpha_v_tangent(ii));
                            theta_step_expr(1+n_contacts+ii) = (1-alpha_q(ii))*(1-alpha_v_normal(ii))*(1-alpha_v_tangent(ii));
                        else
                            theta_step_expr(ii+1) = (1-alpha_q(ii))*(1-alpha_v_normal(ii));
                        end
                    end
                end
                theta_step_expr(1) = alpha_ode;
            else
                % lift and have bilinear terms
                if nonsmooth_switching_fun
                    if n_contacts <= 2
                        for ii = 1:n_contacts
                            alpha_ode = alpha_ode*alpha_qv(ii);
                        end
                        theta_step_expr(1) = alpha_ode;
                    else
                        beta_prod_expr(1) = (alpha_qv(1))*(alpha_qv(2));
                        beta_prod_expr_guess(1) = (alpha_qv(1))*(alpha_qv(2));
                        % lifting terms in between
                        for ii = 3:n_contacts-1
                            beta_prod_expr(ii-1) = beta_prod(ii-2)*(alpha_qv(ii)); % beta_{i} = beta{i-1}*(prod_term_i+1}
                            beta_prod_expr_guess(ii-1) = beta_prod_expr_guess(ii-2)*(alpha_qv(ii)); % this is to have an expression depending only on alpha for the inital guess eval
                        end
                        theta_step_expr(1)= beta_prod(end)*(alpha_qv(n_contacts)); % last lifting term;
                    end
                    % lifting of aux dyn multiplier expressions
                    for ii = 1:n_contacts
                        if friction_exists
                            theta_step_expr(ii+1) = (1-alpha_qv(ii))*(alpha_v_tangent(ii));
                            theta_step_expr(1+n_contacts+ii) = (1-alpha_qv(ii))*(1-alpha_v_tangent(ii));
                        else
                            theta_step_expr(ii+1)=(1-alpha_qv(ii));
                        end
                    end
                else
                    % with two smooth switching functions
                    beta_bilinear_ode_expr = alpha_q.*alpha_v_normal;
                    if friction_exists
                        beta_bilinear_aux_expr = (1-alpha_q).*(1-alpha_v_normal);
                    end

                    if n_contacts <= 2
                        % here no lifting of product terms
                        for ii = 1:n_contacts
                            alpha_ode = alpha_ode*(alpha_q(ii)+alpha_v_normal(ii)-beta_bilinear_ode(ii));
                        end
                        theta_step_expr(1) = alpha_ode;
                    else
                        % here lifting of product terms
                        g_z_tf_beta_prod  = [beta_prod(1) - (alpha_q(1)+alpha_v_normal(1)-beta_bilinear_ode(1))*(alpha_q(2)+alpha_v_normal(2)-beta_bilinear_ode(2))]; % first lifting terms
                        % lifting terms in between
                        for ii = 3:n_contacts-1
                            % beta_{i} = beta{i-1}*(prod_term_i+1}
                            beta_prod_expr(ii-1) = beta_prod(ii-2)*(alpha_q(ii)+alpha_v_normal(ii)-beta_bilinear_ode(ii));
                            beta_prod_expr_guess(ii-1) = beta_prod_expr_guess(ii-2)*(alpha_q(ii)+alpha_v_normal(ii)-beta_bilinear_ode(ii));
                        end
                        % last lifting term;
                        theta_step_expr(1) = beta_prod(end)*(alpha_q(n_contacts)+alpha_v_normal(n_contacts)-beta_bilinear_ode(n_contacts));
                    end
                    % lifting of aux dyn multiplier expressions
                    for ii = 1:n_contacts
                        if friction_exists
                            theta_step_expr(ii+1) =  beta_bilinear_aux(ii)*(alpha_v_tangent(ii));
                            theta_step_expr(1+n_contacts+ii) = beta_bilinear_aux(ii)*(1-alpha_v_tangent(ii));
                        else
                            theta_step_expr(ii+1) = (1-alpha_q(ii))*(1-alpha_v_normal(ii));
                        end
                    end
                end
            end
            % equality constraints in DCS
            g_lift_theta_step = theta_step-theta_step_expr;
            g_lift_beta= beta-[beta_bilinear_ode_expr;beta_bilinear_aux_expr;beta_prod_expr];
            % auxiliary functions to get inital guess for new algebraic variables theta and beta
            g_lift_theta_step_fun  = Function('g_lift_theta_step_fun',{alpha,beta},{theta_step_expr});
            g_lift_beta_fun = Function('g_lift_beta_fun',{alpha},{[beta_bilinear_ode_expr;beta_bilinear_aux_expr;beta_prod_expr_guess]});
            theta_step_all{1} = theta_step;
        end
        n_beta = length(beta);
        n_theta_step = length(theta_step);
        n_z_all = n_z_all + n_beta+n_theta_step;
    case 'CLS'
        n_z_all = n_contacts+n_tangents;
        lambda_normal = define_casadi_symbolic(casadi_symbolic_mode,'lambda_normal',n_contacts);
        y_gap = define_casadi_symbolic(casadi_symbolic_mode,'y_gap',n_contacts);
        g_lift_gap = f_c-y_gap; % lift gap functions f_c(q) = y;
        % Variables for impulse equations
        Lambda_normal = define_casadi_symbolic(casadi_symbolic_mode,'Lambda_normal',n_contacts);
        Y_gap = define_casadi_symbolic(casadi_symbolic_mode,'Y_gap',n_contacts);
        P_vn = define_casadi_symbolic(casadi_symbolic_mode,'P_vn',n_contacts); % pos part of state jump law
        N_vn = define_casadi_symbolic(casadi_symbolic_mode,'N_vn',n_contacts); % neg part of state jump law
        if friction_exists
            % tangetial contact froce (firction force)
            lambda_tangent = define_casadi_symbolic(casadi_symbolic_mode,'lambda_tangent',n_tangents);
            % Impulse varaibles
            Lambda_tangent = define_casadi_symbolic(casadi_symbolic_mode,'Lambda_tangent',n_tangents);
            if isequal(friction_model,'polyhedral')
                gamma_d = define_casadi_symbolic(casadi_symbolic_mode,'gamma_d',n_contacts);
                beta_d = define_casadi_symbolic(casadi_symbolic_mode,'beta_d',n_contacts); % lift friction cone bound
                delta_d = define_casadi_symbolic(casadi_symbolic_mode,'delta_d',n_tangents); % lift lagrangian
                % Impulse varaibles
                Gamma_d = define_casadi_symbolic(casadi_symbolic_mode,'Gamma_d',n_contacts);
                Beta_d = define_casadi_symbolic(casadi_symbolic_mode,'Beta_d',n_contacts); % lift friction cone bound
                Delta_d = define_casadi_symbolic(casadi_symbolic_mode,'Delta_d',n_tangents); % lift lagrangian
                for ii = 1:n_contacts
                    ind_temp = n_t*ii-(n_t-1):n_t*ii;
                    g_lift_friction1 = [g_lift_friction1; beta_d(ii)-(mu(ii)*lambda_normal(ii)- sum(lambda_tangent(ind_temp)))];
                    g_lift_friction2 = [g_lift_friction2; delta_d(ind_temp) - (D_tangent(:,ind_temp)'*v+gamma_d(ii))];
                end
            end
            if isequal(friction_model,'conic')
                gamma = define_casadi_symbolic(casadi_symbolic_mode,'gamma',n_contacts);
                beta = define_casadi_symbolic(casadi_symbolic_mode,'beta',n_contacts);
                % Impulse variables;
                Gamma = define_casadi_symbolic(casadi_symbolic_mode,'Gamma',n_contacts);
                Beta = define_casadi_symbolic(casadi_symbolic_mode,'Beta',n_contacts);
                for ii = 1:n_contacts
                    ind_temp = n_t*ii-(n_t-1):n_t*ii;
                    g_lift_friction1 = [g_lift_friction1; beta(ii)-((mu(ii)*lambda_normal(ii))^2- norm(lambda_tangent(ind_temp))^2)];
                end
                switch conic_friction_switch_detection_mode
                    case 'none'
                        % no extra constraints
                    case 'abs'
                        p_vt = define_casadi_symbolic(casadi_symbolic_mode,'p_vt',n_tangents); % positive parts of tagnetial velocity (for switch detection)
                        n_vt = define_casadi_symbolic(casadi_symbolic_mode,'n_vt',n_tangents); % negative parts of tagnetial velocity (for switch detection)
                        % Impulse
                        P_vt = define_casadi_symbolic(casadi_symbolic_mode,'P_vt',n_tangents); 
                        N_vt = define_casadi_symbolic(casadi_symbolic_mode,'N_vt',n_tangents); 
                    case 'lp'
                        p_vt = define_casadi_symbolic(casadi_symbolic_mode,'p_vt',n_tangents); % positive parts of tagnetial velocity (for switch detection)
                        n_vt = define_casadi_symbolic(casadi_symbolic_mode,'n_vt',n_tangents); % negative parts of tagnetial velocity (for switch detection)
                        alpha_vt = define_casadi_symbolic(casadi_symbolic_mode,'alpha_vt ',n_tangents); % step function of tangential velocities
                        % impulse
                        P_vt = define_casadi_symbolic(casadi_symbolic_mode,'P_vt',n_tangents); 
                        N_vt = define_casadi_symbolic(casadi_symbolic_mode,'N_vt',n_tangents); 
                        Alpha_vt = define_casadi_symbolic(casadi_symbolic_mode,'Alpha_vt ',n_tangents); 
                end
            end
        end
end
g_lift = [g_lift_theta_step;g_lift_beta;g_lift_gap;g_lift_friction1;g_lift_friction2];

%% Collect algebaric varaibles for the specific DCS mode, define initial guess and bounds
% TODO: @Anton: Do the bounds and guess specification already while defining the varaibles?

z_impulse = []; % only for dcs_mode = cls, note that they are evaluated only at left boundary point of at FE
lbz_impulse = [];
ubz_impulse = [];
z_impulse0 = [];
switch dcs_mode
    case 'Stewart'
        % symbolic variables z = [theta;lambda;mu];
        z_all = [vertcat(theta_all{:});vertcat(lambda_all{:});vertcat(mu_all{:})];
        z_switching = [vertcat(lambda_all{:});vertcat(mu_all{:})];
        lbz_all = [0*ones(n_theta,1);0*ones(n_theta,1);-inf*ones(n_sys,1)];
        ubz_all = [inf*ones(n_theta,1);inf*ones(n_theta,1);inf*ones(n_sys,1)];
        % initial guess for z; % solve LP for guess;
        if lp_initialization
            [theta_guess,lambda_guess,mu_guess] = create_lp_based_guess(model);
        else
            theta_guess = initial_theta*ones(n_theta,1);
            lambda_guess = initial_lambda*ones(n_theta,1);
            mu_guess = initial_mu*ones(n_sys,1);
        end
        z0_all = [theta_guess;lambda_guess;mu_guess];
        n_lift_eq = n_sys;
    case 'Step'
        z_all = [alpha;lambda_n;lambda_p;beta;theta_step];
        z_switching = [lambda_n;lambda_p];
        lbz_all = [0*ones(n_alpha,1);0*ones(n_alpha,1);0*ones(n_alpha,1);-inf*ones(n_beta,1);-inf*ones(n_theta_step,1)];
        ubz_all = [ones(n_alpha,1);inf*ones(n_alpha,1);inf*ones(n_alpha,1);inf*ones(n_beta,1);inf*ones(n_theta_step,1)];

        alpha_guess = initial_alpha*ones(n_alpha,1);
        lambda_0_guess = initial_lambda_0*ones(n_alpha,1);
        lambda_1_guess = initial_lambda_1*ones(n_alpha,1);
        beta_guess = initial_beta*ones(n_beta,1);
        theta_step_guess = initial_theta_step*ones(n_theta_step,1);
        % TODO beta guess should exist once custom lifting is implemented
        %                 if pss_lift_step_functions && ~settings.general_inclusion
        if pss_lift_step_functions && ~settings.general_inclusion
            beta_guess = full(g_lift_beta_fun(alpha_guess));
            theta_step_guess = full(g_lift_theta_step_fun(alpha_guess,beta_guess));
        end
        % eval functios for theta_step and beta?
        z0_all = [alpha_guess;lambda_0_guess;lambda_1_guess;beta_guess;theta_step_guess];
        n_lift_eq =length(g_lift);
    case 'CLS'
        z_all = [lambda_normal;y_gap];
        lbz_all = [0*ones(n_contacts,1);0*ones(n_contacts,1)];
        ubz_all = [inf*ones(n_contacts,1);inf*ones(n_contacts,1)];
        z0_all = [ones(n_contacts,1);ones(n_contacts,1)];
        % Impulse
        z_impulse = [Lambda_normal;Y_gap;P_vn;N_vn]; 
        lbz_impulse = [0*ones(n_contacts,1);0*ones(n_contacts,1);0*ones(n_contacts,1);0*ones(n_contacts,1)];
        ubz_impulse = [inf*ones(n_contacts,1);inf*ones(n_contacts,1);inf*ones(n_contacts,1);inf*ones(n_contacts,1)];
        z_impulse0 = [ones(n_contacts,1);ones(n_contacts,1);ones(n_contacts,1);ones(n_contacts,1)];
        if friction_exists
            % tangetial contact froce (firction force)
            z_all = [z_all;lambda_tangent];
            z0_all = [z0_all; ones(n_tangents,1)];
            % Impulse
            z_impulse = [z_impulse;Lambda_tangent]; % only for dcs_mode = cls, note that they are evaluated only at left boundary point of at FE
            z_impulse0 = [z_impulse0; ones(n_tangents,1)];
            % friction aux multipliers
            if isequal(friction_model,'polyhedral')
                % bounds on friction force
                lbz_all = [lbz_all;0*ones(n_tangents,1)];
                ubz_all = [ubz_all;inf*ones(n_tangents,1)];
                % Impulse
                lbz_impulse = [lbz_impulse;0*ones(n_tangents,1)];
                ubz_impulse = [ubz_impulse; inf*ones(n_tangents,1)];

                % polyhedral friction model algebaric variables
                z_all = [z_all;gamma_d;beta_d;delta_d];
                lbz_all = [lbz_all;0*ones(n_contacts,1);0*ones(n_contacts,1);0*ones(n_tangents,1)];
                ubz_all = [ubz_all;inf*ones(n_contacts,1);inf*ones(n_contacts,1);inf*ones(n_tangents,1)];
                z0_all = [z0_all;ones(n_contacts,1);ones(n_contacts,1);ones(n_tangents,1)];
                % Polyhedral friction - collect impulse variables
                z_impulse = [z_impulse;Gamma_d;Beta_d;Delta_d]; 
                lbz_impulse = [lbz_impulse;0*ones(n_contacts,1);0*ones(n_contacts,1);0*ones(n_tangents,1)];
                ubz_impulse = [ubz_impulse ;inf*ones(n_contacts,1);inf*ones(n_contacts,1);inf*ones(n_tangents,1)];
                z_impulse0 = [z_impulse0;ones(n_contacts,1);ones(n_contacts,1);ones(n_tangents,1)];
            end
            if isequal(friction_model,'conic')
                % bounds for friction froce
                lbz_all = [lbz;-inf*ones(n_tangents,1)];
                ubz_all = [ubz;inf*ones(n_tangents,1)];
                % Impulse
                lbz_impulse = [lbz_impulse;-inf*ones(n_tangents,1)];
                ubz_impulse = [ubz_impulse; inf*ones(n_tangents,1)];
                % conic friction model algebaric variables
                z_all = [z_all;gamma;beta];
                lbz_all = [lbz_all;0*ones(n_contacts,1);0*ones(n_contacts,1)];
                ubz_all = [ubz_all;inf*ones(n_contacts,1);inf*ones(n_contacts,1)];
                z0_all = [z0_all; ones(n_contacts,1);ones(n_contacts,1)];
                % Conic impulse
                z_impulse = [z_impulse;Gamma;Beta]; 
                lbz_impulse = [lbz_impulse;0*ones(n_contacts,1);0*ones(n_contacts,1)];
                ubz_impulse = [ubz_impulse;inf*ones(n_contacts,1);inf*ones(n_contacts,1)];
                z_impulse0 = [z_impulse0;ones(n_contacts,1);ones(n_contacts,1)];
                switch conic_friction_switch_detection_mode
                    case 'none'
                        % no extra constraints
                    case 'abs'
                        z_all = [z_all;p_vt;n_vt];
                        lbz_all = [lbz_all;0*ones(2*n_tangents,1)];
                        ubz_all = [ubz_all;inf*ones(2*n_tangents,1)];
                        z0_all = [z0_all; ones(2*n_tangents,1)];
                        % Impulse
                        z_impulse = [z_impulse;P_vt;N_vt]; 
                        lbz_impulse = [lbz_impulse;0*ones(2*n_tangents,1)];
                        ubz_impulse = [ubz_impulse;inf*ones(2*n_tangents,1)];
                        z_impulse0 = [z_impulse0; ones(2*n_tangents,1)];
                    case 'lp'
                        z_all = [z_all;p_vt;n_vt;alpha_vt];
                        lbz_all = [lbz_all;0*ones(3*n_tangents,1)];
                        ubz_all = [ubz_all;inf*ones(2*n_tangents,1);1*ones(n_tangents,1)];
                        z0_all = [z0_all; 1*ones(2*n_tangents,1);0.5*ones(n_tangents,1)];
                        % Impulse
                        z_impulse = [z_impulse;P_vt;N_vt;Alpha_vt]; 
                        lbz_impulse = [lbz_impulse;0*ones(3*n_tangents,1)];
                        ubz_impulse = [ubz_impulse;inf*ones(2*n_tangents,1);*ones(n_tangents,1)];
                        z_impulse0 = [z_impulse0; 1*ones(2*n_tangents,1);0.5*ones(n_tangents,1)];
                end
            end
        end
        n_z_all = length(z_all);
end

%% Add user provided algebraic
z_all = vertcat(z_all,z);
z0_all = [z0_all;z0];
lbz_all = [lbz_all;lbz];
ubz_all = [ubz_all;ubz];
n_z_all = n_z_all + n_z;
%% TODO %%%%%%%%%%%%%%
%%%%%% HERE ENDS SECOND FUNCTION THAT DEFINES THE MODEL VARIABLES. THE USER
%%%%%% MIGHT USE THEM IN THE CONSTRAINTS (BUT DO WE GO THROUGH ANOTHER CHECK OF G_INEQ AND CO?
% ----> THE NEXT PART IS THE THIRD AND LAST FUNCTION THAT DEFINES THE CASADI FUNCTIONS AND DIFFERENTIAL EQUATIONS

%% Model functions of the DCS mode
% if f_x doesnt exist we generate it from F
% if it does we are in expert mode. TODO name.
% Define differential equations
if ~isfield(model, 'f_x')
    f_x = zeros(n_x,1);
    % rhs of ODE;
    for ii = 1:n_sys
        switch dcs_mode
            case 'Stewart'
                f_x = f_x + F{ii}*theta_all{ii};
            case 'Step'
                f_x = f_x + F{ii}*theta_step_all{ii};
            case 'DCS'
                if ~lift_velocity_state
                    if friction_exists
                        switch friction_model
                            case 'Conic'
                                    F_v = inv(M)*(f_v+J_normal*lambda_normal+J_tangent*lambda_tangent);
                            case 'Polyhedral'
                                F_v = inv(M)*(f_v+J_normal*lambda_normal + D_tangent*lambda_tangent);
                        end
                    else
                        F_v = inv(M)*(f_v+J_normal*lambda_n);
                    end
                     f_x = [v;F_v];
                else
                    z_v = define_casadi_symbolic(casadi_symbolic_mode,['z_v'],n_q);
                    z_all = [z_all;z_v];
                    n_z_all = n_z_all+n_q;
                    lbz_all = [lbz_all;-inf*ones(n_q,1)];
                    ubz_all = [ubzz_all;inf*ones(n_q,1)];
                    z0_all = [z0_all;x0(n_q+1:end)];

                    f_x = [v;z_v];
                    if friction_exists
                        switch friction_model
                            case 'Conic'
                                g_lift_v = M*z_v -(f_v +J_normal*lambda_normal + J_tangent*lambda_tangent);
                            case 'Polyhedral'
                                g_lift_v =  M*z_v -(f_v + J_normal*lambda_normal + D_tangent*lambda_tangent);
                        end
                                                
                    else
                        g_lift_v =  M*z_v -(f_v+J_normal*lambda_n);
                    end
                    g_lift = [g_lift;g_lift_v];
                end
        end
    end
end

% Define algebraic equations

g_switching = []; % collects switching function algebraic equations 0 = g_i(x) - \lambda_i - e \mu_i
g_convex = []; % equation for the convex multiplers 1 = e' \theta
f_comp_residual = 0; % the orthogonality conditions diag(\theta) \lambda = 0.
lambda00_expr =[];
for ii = 1:n_sys
    switch dcs_mode
        case 'Stewart'
            % basic algebraic equations and complementarity condtions of the DCS
            % (Note that the cross complementarities are later defined when the discrete
            % time variables for every IRK stage in the create_nlp_nosnoc function are defined.)
            % g_ind_i - lambda_i + mu_i e_i = 0; for all i = 1,..., n_sys
            % lambda_i'*theta_i = 0; for all i = 1,..., n_sys
            % lambda_i >= 0;    for all i = 1,..., n_sys
            % theta_i >= 0;     for all i = 1,..., n_sys
            % Gradient of Lagrange Function of indicator LP
            g_switching = [g_switching; g_Stewart{ii}-lambda_all{ii}+mu_all{ii}*e_ones_all{ii}];
            g_convex = [g_convex;e_ones_all{ii}'*theta_all{ii}-1];
            lambda00_expr = [lambda00_expr; g_Stewart{ii}- min(g_Stewart{ii})];
            f_comp_residual = f_comp_residual + lambda_all{ii}'*theta_all{ii};
        case 'Step'
            % c_i(x) - (lambda_1_i-lambda_0_i)  = 0; for all i = 1,..., n_sys
            % lambda_0_i'*alpha_i  = 0; for all i = 1,..., n_sys
            % lambda_1_i'*(e-alpha_i)  = 0; for all i = 1,..., n_sys
            % lambda_0_i >= 0;    for all i = 1,..., n_sys
            % lambda_1_i >= 0;    for all i = 1,..., n_sys
            % alpha_i >= 0;     for all i = 1,..., n_sys
            g_switching = [g_switching;c{ii}-lambda_1_all{ii}+lambda_0_all{ii}];
            f_comp_residual = f_comp_residual + lambda_0_all{ii}'*alpha_all{ii}+lambda_1_all{ii}'*(ones(n_c_sys(ii),1)-alpha_all{ii});
            %             lambda00_expr = [lambda00_expr; max(c{ii},0); -min(c{ii}, 0)];
            lambda00_expr = [lambda00_expr; -min(c{ii}, 0); max(c{ii},0)];
        case 'DCS'
            % Algebraic equations
%                 g_z_comp1 = [];
%                 g_z_comp2 = [];
%                 if friction_exists
%                     if isequal(friction_model,'Polyhedral')
%                         for ii = 1:n_contacts
%                             ind_temp = n_t*ii-(n_t-1):n_t*ii;
%                             g_z_comp1 = [lambda_t(ind_temp);gamma_d(ii)];
%                             g_z_comp2 = [D_t(:,ind_temp)'*v-gamma_d(ii);mu(ii)*lambda_n(ii)-sum(lambda_t(ind_temp))];
%                         end
%                     end
%                     if isequal(friction_model,'Conic')
%                         for ii = 1:n_contacts
%                             ind_temp = n_t*ii-(n_t-1):n_t*ii;
%                             g_z = [g_z;...
%                                 beta(ii)-( (mu(ii)*lambda_n(ii))^2-norm(lambda_t(ind_temp))^2);
%                                 -J_t(:,ind_temp)'*v-2*gamma(ii)*(lambda_t(ind_temp)+kappa_reg_tan*f_c(ii))];
%                             g_z_comp1 = [g_z_comp1;gamma(ii)];
%                             g_z_comp2 = [g_z_comp2;beta(ii)];
% 
%                             if isequal(conic_friction_switch_detection_mode,'abs')|| isequal(conic_friction_switch_detection_mode,'lp')
%                                 % The equality constraint is the same in both variantes
%                                 g_z = [g_z;J_t(:,(ind_temp))'*v-(p_vt(ind_temp)-n_vt(ind_temp))];
%                             end
%                         end
%                     end
%                 end
%                 lbg_z = zeros(length(g_z),1);
%                 ubg_z = zeros(length(g_z),1);
%                 ubg_z(1:n_contacts) = inf;
 % contact normal
%                 invM = inv(M);
%                 M_fun = Function('M_fun', {x}, {M});
%                 invM_fun = Function('invM_fun', {x}, {invM});
%                 f_c_fun = Function('f_c_fun', {x}, {f_c});
%                 J_n_fun = Function('J_n_fun', {x}, {J_n});
%                 if friction_exists
%                     if isequal(friction_model,'Conic')
%                         J_t_fun = Function('J_t_fun', {x}, {J_t});
%                     else
%                         D_t_fun = Function('D_t_fun', {x}, {D_t});
%                     end
%                 end

    end
end

%% Lifting of forces in rigid bodies (either with time freezing or the dcs system)
% the r.h.s of M(q)\ddot{q} = F(q,\dot{q},u) into  M{q}z_v-f(q,dor{q},u)=0; \ddot{q} = z_v
% TODO: rename time_freezing_lift_forces into lift_velocity_state
g_lift_forces = [];
if time_freezing && time_freezing_lift_forces % TODO Is this broken with parameters/v_global
    f_v = f_x(n_q+1:2*n_q);
    if n_u > 0
        f_v_fun = Function('f_v_fun',{x,z_all,u,v_global},{f_v});
        z0_forces = full(f_v_fun(x0,z0_all,u0,v0_global));
    else % TODO remove this?
        f_v_fun = Function('f_v_fun',{x,z_all,v_global},{f_v});
        z0_forces = full(f_v_fun(x0,z0_all,v0_global));
    end
    z_forces = define_casadi_symbolic(casadi_symbolic_mode,'z_forces',n_q);
    z_all = [z_all;z_forces];
    n_z_all = n_z_all+n_q;
    z0_all = [z0_all;z0_forces];
    lbz_all = [lbz_all;-inf*ones(n_q,1)];
    ubz_all = [ubz_all;inf*ones(n_q,1)];
    g_lift_forces = [M*z_forces - f_v]; % lifting function
    % new simple dynamics after lifting
    f_x = [f_x(1:n_q); z_forces; f_x(end-n_quad:end)];
end

if lift_velocity_state

end

%%  collect all algebraic equations
g_lp = [g_switching;g_convex;g_lift];
g_z_all = [g_lp;g_z;g_lift_forces];
n_algebraic_constraints = length(g_z_all);

%% CasADi functions for indicator and region constraint functions
% model equations
% TODO: @ Anton make this to be a function of v_global as well (and test on telescop arm example)
g_Stewart_fun = Function('g_Stewart_fun',{x,p},{g_ind_vec});
c_fun = Function('c_fun',{x,p},{c_all});
dot_c = c_all.jacobian(x)*f_x;

f_x_fun = Function('f_x_fun',{x,z_all,u,p,v_global},{f_x,f_q});
g_z_all_fun = Function('g_z_all_fun',{x,z_all,u,p,v_global},{g_z_all}); % lp kkt conditions without bilinear complementarity term (it is treated with the other c.c. conditions)

g_switching_fun = Function('g_switching_fun', {x,z_switching,u,p}, {g_switching});
dot_c_fun = Function('c_fun',{x,z_all,u,p},{dot_c}); % total time derivative of switching functions

model.lambda00_fun = Function('lambda00_fun',{x,p_global},{lambda00_expr});

J_cc_fun = Function('J_cc_fun',{z_all},{f_comp_residual});
f_q_T_fun = Function('f_q_T',{x,p,v_global},{f_q_T});

%%  CasADi functions for lest-square objective function terms
f_lsq_x_fun = Function('f_lsq_x_fun',{x,x_ref,p},{f_lsq_x});
if n_u > 0
    f_lsq_u_fun = Function('f_lsq_u_fun',{u,u_ref,p},{f_lsq_u});
end
f_lsq_T_fun = Function('f_lsq_T_fun',{x,x_ref_end,p_global},{f_lsq_T});

%% Initial guess for state derivatives at stage points
if isequal(irk_representation,'differential')
    if simple_v0_guess
        v0 = zeros(n_x,1);
    else
        [v0,~] = (f_x_fun(x0,z0_all,u0,[p_global_val; p_time_var_val(:,1)], v0_global));
        v0 = full(v0);
    end
    model.v0 = v0;
end
%% Collect Outputs
model.lbx = lbx;
model.ubx = ubx;
model.lbu = lbu;
model.ubu = ubu;
if n_u > 0
    model.u0 = u0;
end
model.z0_all = z0_all;
model.lbz_all = lbz_all;
model.ubz_all = ubz_all;

model.z0 = z0;
model.lbz = lbz;
model.ubz = ubz;

model.v_global = v_global;
model.v0_global = v0_global;
model.lbv_global = lbv_global;
model.ubv_global = ubv_global;

model.g_path_constraint = g_path_constraint;
if g_path_constraint
    model.g_path_lb = g_path_lb;
    model.g_path_ub = g_path_ub;
    model.g_path_fun = g_path_fun;
end

model.g_comp_path_constraint = g_comp_path_constraint;
if g_comp_path_constraint
    model.g_comp_path_fun = g_comp_path_fun;
end

model.terminal_constraint = terminal_constraint;
if terminal_constraint
    model.g_terminal_lb = g_terminal_lb;
    model.g_terminal_ub = g_terminal_ub;
    model.g_terminal_fun = g_terminal_fun;
end

if exist('g_lift_theta_step_fun')
    model.g_lift_theta_step_fun = g_lift_theta_step_fun;
end
if exist('g_lift_beta_fun')
    model.g_lift_beta_fun = g_lift_beta_fun;
end

% CasADi Expressions
model.f_x = f_x;
model.f_q = f_q;
model.g_switching = g_switching;
model.g_z_all = g_z_all;
model.f_q_T = f_q_T;
% CasADi Functions
model.f_x_fun = f_x_fun;
model.g_z_all_fun = g_z_all_fun;
model.g_switching_fun = g_switching_fun;
model.f_q_T_fun = f_q_T_fun;

model.J_cc_fun = J_cc_fun;
model.g_Stewart_fun = g_Stewart_fun;
model.c_fun = c_fun;
model.dot_c_fun = dot_c_fun;
%
% % Model Dimensions;
model.n_x = n_x;
model.n_z_all = n_z_all;
model.n_u = n_u;
model.n_sys = n_sys;

model.e_alpha = e_alpha;

model.m_vec = m_vec;
model.m_ind_vec = m_ind_vec;
model.n_theta = n_theta;
model.n_lambda = n_lambda;
model.n_algebraic_constraints = n_algebraic_constraints;
model.n_lift_eq  = n_lift_eq;

model.n_c_sys = n_c_sys;
model.n_alpha = n_alpha;
model.n_beta = n_beta;
model.n_theta_step = n_theta_step;
model.n_lambda_0 = n_lambda_n;
model.n_lambda_1 = n_lambda_p;

% least square functions and references
model.f_lsq_x_fun = f_lsq_x_fun;
model.x_ref_val = x_ref_val;
if n_u>0
    model.f_lsq_u_fun = f_lsq_u_fun;
end
model.u_ref_val = u_ref_val;
model.f_lsq_T_fun = f_lsq_T_fun;
model.x_ref_end_val = x_ref_end_val;

% global parameters
model.p_global = p_global;
model.p_global_val = p_global_val;

% time varying parameters
% TODO maybe make these functions and actually optimization vars. (actually this might just be algebraic vars)
model.p_time_var = p_time_var;
model.p_time_var_stages = p_time_var_stages;
model.p_time_var_val = p_time_var_val;

model.p_dyn = [p_global, p_time_var_stages];
%% collect all dimensions in one sperate struct as it is needed by several other functions later.
dimensions.N_stages = N_stages;
dimensions.N_finite_elements = N_finite_elements;
dimensions.n_x = n_x;
dimensions.n_f = n_f;
dimensions.n_u = n_u;
dimensions.n_z_all = n_z_all;
dimensions.n_z = n_z;
dimensions.n_s = n_s;
dimensions.n_theta = n_theta;
dimensions.n_sys = n_sys;
dimensions.m_vec = m_vec;
dimensions.m_ind_vec = m_ind_vec;
dimensions.n_c_sys = n_c_sys;
dimensions.n_alpha = n_alpha;
dimensions.n_beta = n_beta;
dimensions.n_theta_step = n_theta_step;
dimensions.n_lambda_0 = n_lambda_n;
dimensions.n_lambda_1 = n_lambda_p;
dimensions.n_f_sys = n_f_sys;
dimensions.n_p_global = n_p_global;
dimensions.n_p_time_var = n_p_time_var;
model.dimensions = dimensions;

end


