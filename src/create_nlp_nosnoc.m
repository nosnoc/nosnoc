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
function [solver, solver_initialization, model, settings] = create_nlp_nosnoc(varargin)
% This functions creates the solver instance for the OCP discretized with FESD (or time-stepping IRK scheme).
% The discretization results in an MPCC which can be solved by various
% reformulations, see below.
% -------------------

% TODO: split problem and solver creation!

%% Import CasADi in the workspace of this function
import casadi.*
%% Read data
model = varargin{1};
settings = varargin{2};
%% Reformulation of the PSS into a DCS
[model,settings] = model_reformulation_nosnoc(model,settings);

settings.create_butcher_tableu(model);

%% Formulate the NLP / Main Discretization loop
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
model.g = g;
model.w = w;
model.p = p;
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
if settings.print_level > 5
    problem.print();
end

%% Model update: all index sets and dimensions
% TODO: Maybe just return the problem, currently trying not to break compatibility for now.
model.ind_x = [problem.ind_x0.'; flatten_ind(problem.ind_x)];
model.ind_v = sort(flatten_ind(problem.ind_v));
model.ind_z_all = problem.ind_z_all; %TODO fix this by breaking compat
model.ind_u = problem.ind_u;
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

end
