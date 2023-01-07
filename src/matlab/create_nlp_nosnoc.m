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
% function [solver,solver_initialization, model,settings] = create_nlp_nosnoc(model,settings)
function [varargout] = create_nlp_nosnoc(varargin)
% This functions creates the solver instance for the OCP discretized with FESD (or time-stepping IRK scheme).
% The discretization results in an MPCC which can be solved by various
% reformulations, see below.
% -------------------
% Brief MPCC Wiki
% There are several possible MPCC Solution strategies avilable, by setting mpcc_mode to :
% 'direct' - treat complementarity conditions directly in the NLP, the bilinear term is tread as an inequality constraint.
% 'Scholtes_eq' - Smooth the complementarity conditions, Scholtes' smoothing.
% 'Scholtes_ineq' - Relax the complementarity conditions, Scholtes' relaxation.
% 'ell_1_penalty' - \ell_1 penalty, penalize the sum of all bilinear terms in the objective
% 'elastic_ineq' - \ell_infty elastic mode, upper bound all bilinear term with a positive slack, and penalize the slack in the objective.
% 'elastic_eq' - \ell_infty elastic mode, equate all bilinear term to a positive slack, and penalize the slack in the objective.
% 'elastic_two_sided' - \ell_infty, same as 'elastic_ineq' but two sided.
% 'elastic_ell_1_ineq' - \ell_1, elastic mode but penalize ell_1 norm of complementarities
% 'elastic_ell_1_eq' - \ell_1, elastic mode but penalize ell_1 norm of complementarities
% 'elastic_ell_1_two_sided' - \ell_1, elastic mode but penalize ell_1 norm of complementarities

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

%%  Butcher Tableu
% TODO clean this up.
switch settings.irk_representation
  case 'integral'
    [B,C,D,tau_root] = generate_butcher_tableu_integral(model.dimensions.n_s,settings.irk_scheme);
    if tau_root(end) == 1
        right_boundary_point_explicit  = 1;
    else
        right_boundary_point_explicit  = 0;
    end
    settings.B_irk = B;
    settings.C_irk = C;
    settings.D_irk = D;
  case {'differential', 'differential_lift_x'}
    [A_irk,b_irk,c_irk,order_irk] = generate_butcher_tableu(model.dimensions.n_s,settings.irk_scheme);
    if c_irk(end) <= 1+1e-9 && c_irk(end) >= 1-1e-9
        right_boundary_point_explicit  = 1;
    else
        right_boundary_point_explicit  = 0;
    end
    settings.A_irk = A_irk;
    settings.b_irk = b_irk;
  otherwise
    error('Choose irk_representation either: ''integral'' or ''differential''')
end
settings.right_boundary_point_explicit = right_boundary_point_explicit;


%% Formulate the NLP / Main Discretization loop
% TODO cleanup steps:
%      - Create primal variables all at once.
%      - Separate sections into separate functions operating on the `problem` struct/class
%      - time variables should probably not just be lumped into the state, for readability.
%      - remove index in symbolic variable defintions and add instructive
%        names, e.g., Uk -> U,  h_ki -> h_fe, X_ki_stages ->  X_rk_stages
%      - provide instructive names for terminal constraint relaxations
%      - provide more instructive names for cross_comp (match python)
problem = NosnocProblem(settings, model.dimensions, model);

%% CasADi Functions for objective complementarity residual
w = problem.w; % vectorize all variables, TODO: again, further cleanup necessary
g = problem.g; % vectorize all constraint functions
p = problem.p;
J_fun = problem.cost_fun;
comp_res = problem.comp_res;
comp_res_fesd = problem.comp_fesd;
comp_res_std = problem.comp_std;

%% NLP Solver
prob = struct('f', problem.cost, 'x', w, 'g', g,'p',p);
solver = nlpsol(settings.solver_name, 'ipopt', prob, settings.opts_ipopt);

%% Define CasADi function for the switch indicator function.
nu_fun = Function('nu_fun', {w,p},{problem.nu_vector});

%% Outputs
model.prob = prob;
model.problem = problem;
model.solver = solver;
model.g =  g;
model.w =  w;
model.p =  p;
model.J = problem.cost;
model.J_fun = J_fun;
model.comp_res = comp_res;
model.comp_res_fesd = comp_res_fesd;
model.comp_res_std = comp_res_std;
model.nu_fun = nu_fun;

% create CasADi function for objective gradient.
nabla_J = problem.cost.jacobian(model.w);
nabla_J_fun = Function('nabla_J_fun', {w,p},{nabla_J});
model.nabla_J = nabla_J;
model.nabla_J_fun = nabla_J_fun;

% TODO: make member function
if settings.print_level > 1
    problem.print();
end

%% Model update: all index sets and dimensions
% TODO: Maybe just return the problem, currently trying not to break compatibility for now.
model.ind_x = [problem.ind_x0.'; flatten_ind(problem.ind_x)];
model.ind_v = flatten_ind(problem.ind_v);
model.ind_z = problem.ind_z; %TODO fix this by breaking compat
model.ind_u = flatten_ind(problem.ind_u);
model.ind_h = flatten_ind(problem.ind_h);
model.ind_sot = flatten_ind(problem.ind_sot);
model.ind_t_final  = problem.ind_t_final;
model.p_val = problem.p0;

%% Store solver initialization data
solver_initialization.w0 = problem.w0;
solver_initialization.lbw = problem.lbw;
solver_initialization.ubw = problem.ubw;
solver_initialization.lbg = problem.lbg;
solver_initialization.ubg = problem.ubg;

%% Output
varargout{1} = solver;
varargout{2} = solver_initialization;
varargout{3} = model;
varargout{4} = settings;
end
