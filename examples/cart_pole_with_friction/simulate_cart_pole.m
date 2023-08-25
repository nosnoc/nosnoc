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

% simple simulation
F_friction = 2.0;
model = get_cart_pole_with_friction_model(0, F_friction);

Nsim = 30;
Tsim = 5;
h = Tsim/Nsim;

%% create casadi integrator
import casadi.*
p_integrator = vertcat(model.u, model.p);
ode = struct('x', model.x, 'p', p_integrator, 'ode', model.f_expl_ode);
Phi = integrator('F', 'idas', ode, struct('tf', h));

%% simulation loop
nx = length(model.x);
p_val = ones(length(p_integrator), 1);
p_traj = repmat(p_val, 1, Nsim);
p_traj(1, Nsim/2:end) = -2;
xcurrent = zeros(nx, 1);
simX = zeros(nx, Nsim+1);

simX(:, 1) = xcurrent;
for i = 1:Nsim
    out = Phi('x0', xcurrent, 'p', p_traj(:, i));
    xcurrent = full(out.xf);
    simX(:, i+1) = xcurrent;
end

t_grid = linspace(0, Tsim, Nsim+1);
results = struct('x', simX, 't_grid', t_grid, 't_grid_u', t_grid, 'u', p_traj(1, :));
problem_options.h_k = h;
x_ref = zeros(4, 1);
plot_cart_pole_trajectory(results, model, x_ref);
