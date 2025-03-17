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
%% Description
% In this example a swing up of a pendulum on a cart subject to friction is
% treated. 
% The script allows to solve the problem with a smoothed friction model and
% standard direct collocation, to compare to a nonsmooth friction model and
% nosnoc please run cart_pole_nonsmooth 
% This allows to explore the drawbacks of naive smoothing, e.g., if started
% with a very low smoothing parameter.
% For more details see: https://www.syscop.de/files/2023ss/nonsmooth_school/ex1_sol.pdf
close all
clear all
%% Model
F_friction = 2; % Friction force amplitude

% model
model = get_cart_pole_with_friction_model(true, F_friction);
x_ref = [0; 180/180*pi; 0; 0]; % target position

% Discretization options
problem_options = nosnoc.Options();
problem_options.T = 5;  % Time horizon
problem_options.rk_scheme = RKSchemes.RADAU_IIA;
%problem_options.rk_representation = 'differential';
n_s = 3;
problem_options.n_s = n_s;
problem_options.dcs_mode = 'Stewart';
problem_options.N_stages = 2; % number of control intervals
N_fe = 2;
problem_options.N_finite_elements = N_fe; % number of finite element on every control interval
problem_options.cross_comp_mode = 'FE_FE';

% solver options
solver_options = nosnoc.reg_homotopy.Options();
solver_options.N_homotopy = 15;
solver_options.complementarity_tol = 1e-13;
solver_options.sigma_N = 1e-13;
%solver_options.solver = 'fatrop';

%mpecopt options
%solver_options = mpecopt.Options();
%solver_options.rho_TR_phase_i_init = 1;

% other linear solvers require installation, check https://www.hsl.rl.ac.uk/catalogue/ and casadi.org for instructions
% solver_options.opts_casadi_nlp.ipopt.linear_solver = 'ma27';

% create solver and solve problem
ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);
ocp_solver.solve();

close all;
nlp = ocp_solver.discrete_time_problem.solver.nlp;
dyn = [];
alg = [];
int = [];
w = [nlp.w.x(0,0,2);nlp.w.lambda(0,0,2)];
w_alt = [nlp.w.x(0,0,2);nlp.w.lambda(0,0,2)];
for ii=1:N_fe
    x = [];
    lam = [];
    theta = [];
    mu = [];
    for jj=1:n_s
        dyn = [dyn;nlp.g.dynamics(1,ii,jj)];
        alg = [alg;nlp.g.algebraic(1,ii,jj)];
        int = [int;nlp.g.dynamics(1,ii,jj);nlp.g.algebraic(1,ii,jj)];
        w = [w;nlp.w.x(1,ii,jj);nlp.w.lambda(1,ii,jj);nlp.w.theta(1,ii,jj);nlp.w.mu(1,ii,jj)];
        x = [x;nlp.w.x(1,ii,jj)];
        lam = [lam;nlp.w.lambda(1,ii,jj)];
        theta = [theta;nlp.w.theta(1,ii,jj)];
        mu = [mu;nlp.w.mu(1,ii,jj)];
    end
    w = [w;nlp.w.h(1,ii)]
    w_alt = [w_alt;x;lam;theta;mu;nlp.w.h(1,ii)];
end
w = [w; nlp.w.u(1)];;
w_alt = [w_alt;nlp.w.u(1)];
figure
spy(dyn.jacobian(w))
title("dyn,interleaved")
figure
spy(dyn.jacobian(w_alt))
title("dyn,blocks")
figure
spy(alg.jacobian(w))
title("alg,interleaved")
figure
spy(alg.jacobian(w_alt))
title("alg,blocks")
figure
spy(int.jacobian(w))
title("int,interleaved")
figure
spy(int.jacobian(w_alt))
title("int,blocks")
