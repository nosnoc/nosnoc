function [varargout] = time_freezing_kinematics_iteration(varargin)
%% preproces settings
model = varargin{1};
settings = varargin{2};
if nargin > 2
    w0 = varargin{3};
end
import casadi.*

if settings.virtual_forces_convex_combination 
    remove_psi_vf = 0;
else
    remove_psi_vf = 1;
end

if settings.virtual_forces
    remove_u_virtual = 0;
else
    n_u_no_vf = length(model.u);
    remove_u_virtual = 1;
end

%% Updated settings - virtual forces
settings.virtual_forces = 1;
settings.tighthen_virtual_froces_bounds = 0; % squeez bounds for virual forces to zero
settings.penalize_virtual_forces = 1;  % increasing qudratic penalty for virtual forces
settings.virtual_forces_convex_combination = 1;  % 1- convex combination between kinematics and true dynamics, 0 -
% take from the input/default
% settings.virtual_forces_in_every_mode = 1; % 0 -is it just in uncondraind dynamics (nonsmoothnes presevred in convex mode), 1-it is in every pss mode (smooth kinematics ode in convex mode)
% settings.virtual_forces_in_every_mode = 0; (% usually more difficult problem)
settings.virtual_forces_parametric_multipler = 0; % 1- multiplier is external parameter, 0 - optimization variable
settings.M_virtual_forces = 1e2; % bound for virtual forces

n_q = length(model.f);
% remove dynamics
model.f = zeros(n_q,1);
model.M = eye(n_q);
% remove objective
% model.f_q = 0;
% model.f_q_T = 0;
%% Create NLP element
tic
[solver,solver_initalization,model,settings] = create_nlp_nosnoc(model,settings);
solver_generating_time = toc;
if settings.print_level >=2
    fprintf('Kinematics solver generated in in %2.2f s. \n',solver_generating_time);
end

%% Check provided initial guess
if exist("w0",'var')
    if length(w0) == length(solver_initalization.w0)
        solver_initalization.w0 = w0;
    else
        fprintf('Provided user guess does not have the appropiate dimension, it should be a vector of length %d, the provided vectors has a length of %d. \n',length(solver_initalization.w0),length(w0));
        fprintf('Taking the generated default initial guess... \n');
    end
end

%% Update bounds
solver_initalization.lbw(model.ind_vf) = 1;
solver_initalization.ubw(model.ind_vf) = 1;
%% read data
unfold_struct(solver_initalization,'caller');
unfold_struct(model,'caller');
unfold_struct(settings,'caller');

%% solve nlp
tic
results = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg,'p',p_val);
cpu_time_iter = toc;
w0 = full(results.x);
complementarity_iter = full(comp_res(w0));
vf_resiudal = full(J_virtual_froces_fun(w0));
if print_level >= 3
    fprintf('-----------------------------------------------------------------------------------------------\n');
    fprintf('Solved problem with kinematic model.\n');
    fprintf('Complementarity resiudal: %2.2e.\n',complementarity_iter);
    fprintf('CPU time of iteration: %2.2f s.\n',cpu_time_iter);
    fprintf('Virtual forces residual: %2.2e.\n',vf_resiudal);
    fprintf('-----------------------------------------------------------------------------------------------\n');
end
%% post process
w0_unchanged = w0;

if remove_psi_vf 
    w0(ind_vf) = [];
end

if remove_u_virtual
    ind_u = reshape(ind_u,n_u,N_stages);
    ind_u(1:n_u_no_vf,:) = [];
    ind_u = ind_u(:);
    w0(ind_u) = [];
end
%% output
varargout{1} = w0;
varargout{2} = cpu_time_iter;
varargout{3} = w0_unchanged;

end