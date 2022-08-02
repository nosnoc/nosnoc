%    This file is part of NOSNOC.
%
%    NOSNOC -- A software for NOnSmooth Numerical Optimal Control.
%    Copyright (C) 2022 Armin Nurkanovic, Moritz Diehl (ALU Freiburg).
%
%    NOSNOC is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 3 of the License, or (at your option) any later version.
%
%    NOSNOC is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with NOSNOC; if not, write to the Free Software Foundation,
%    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%
function [varargout] = create_nlp_nosnoc_opti(varargin)
% Info on this function.
%% read data
model = varargin{1};
settings = varargin{2};
%% CasADi
import casadi.*
%% Reformulation of the PSS into a DCS
[settings] = refine_user_settings(settings);
[model,settings] = model_reformulation_nosnoc(model,settings);

%% Fillin missing settings with default settings
[settings] = fill_in_missing_settings(settings,model);

%% Load user settings and model details
unfold_struct(settings,'caller')
unfold_struct(model,'caller');
%% Bounds on step-size (if FESD active)

%% Butcher Tableu (differential and integral representation)
switch irk_representation
    case 'integral'
        [B,C,D,tau_root] = generatre_butcher_tableu_integral(n_s,irk_scheme);
        if tau_root(end) == 1
            right_boundary_point_explicit  = 1;
        else
            right_boundary_point_explicit  = 0;
        end
    case 'differential'
        [A_irk,b_irk,c_irk,order_irk] = generatre_butcher_tableu(n_s,irk_scheme);
        if c_irk(end) <= 1+1e-9 && c_irk(end) >= 1-1e-9
            right_boundary_point_explicit  = 1;
        else
            right_boundary_point_explicit  = 0;
        end
    otherwise
        error('Choose irk_representation either: ''integral'' or ''differential''')
end
%% Elastic Mode variables (do this below)

%% Time optimal control (time transformations)

% clc; clear all; close all;
%% Discretization
n_s = 2;
N_stg = 20; % number of control intervals
N_FE = 3;

%% Butcher tables, differential and integral representation
tau = collocation_points(n_s, 'radau');
[C,D,B] = collocation_coeff(tau);

%% Example OCP
% Time horizon
T = 10;
% Declare model variables
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1; x2];
x0 = [0;1];
lbx = [-2;-2];
ubx = [2;2];
n_x = 2;
u = SX.sym('u');
lbu = -1;
ubu = 1;
u0 = 1;
n_u = 1;
% algebraic variables
z = SX.sym('z');
n_z = 1;
z0 = 1;
lbz = -inf;
ubz = inf;


% Model equations
f_x = [(1-x2^2)*x1 - x2 + u; x1-0*tanh(x2/1e-2)];
g_z = x1^2+x2^2-z;
% Objective term
f_q = x1^2 + x2^2 + u^2;
% Continuous time dynamics
f_x_fun = Function('f_x_fun', {x, z, u}, {f_x, f_q});
g_z_fun = Function('g_z_fun', {x, z, u}, {g_z});
% Control discretization
h = T/(N_stg*N_FE);
%% Start with an empty NLP
opti = Opti();
J = 0;
% "Lift" initial conditions
Xk = opti.variable(n_x);
opti.subject_to(Xk==x0);
opti.set_initial(Xk, x0);
% Index sets
ind_total = [1:n_x];
ind_x = [1:n_x];
ind_u = [];
ind_z = [];
% Collect all states/controls
Xs = {Xk};
Xs_extended = {Xk};
Us = {};
Zs = {};

%% Formulate the NLP
for k=0:N_stg-1
    %% New NLP variable for the control
    Uk = opti.variable(n_u);
    Us{end+1} = Uk;
    opti.subject_to(lbu <= Uk <=ubu);
    opti.set_initial(Uk, 0);
    ind_u = [ind_u,ind_total(end)+1:ind_total(end)+n_u];
    ind_total = [ind_total,ind_total(end)+1:ind_total(end)+n_u];
    %% Loop over finite elements for the current control interval
    for i = 0:N_FE-1
        %& Variables at stage points
        Xc = opti.variable(n_x, n_s);
        Xs_extended{end+1} = Xc;
        opti.subject_to(lbx <= Xc <=ubx);
        opti.set_initial(Xc, repmat(x0,1,n_s));
        ind_x = [ind_x,ind_total(end)+1:ind_total(end)+n_x*n_s];
        ind_total = [ind_total,ind_total(end)+1:ind_total(end)+n_x*n_s];

        Zc = opti.variable(n_z, n_s);
        Zs{end+1} = Zc;

        opti.subject_to(lbz <= Zc <=ubz);
        opti.set_initial(Zc, repmat(z0,1,n_s));
        ind_z = [ind_z,ind_total(end)+1:ind_total(end)+n_z*n_s];
        ind_total = [ind_total,ind_total(end)+1:ind_total(end)+n_z*n_s];

        % Evaluate ODE right-hand-side at all stage points
        [f_x, f_q] = f_x_fun(Xc, Zc, Uk);
        [g_z] = g_z_fun(Xc, Zc, Uk);
        % Add contribution to quadrature function
        J = J + f_q*B*h;
        % Get interpolating points of collocation polynomial
        Z = [Xk Xc];
        % Get slope of interpolating polynomial (normalized)
        Pidot = Z*C;
        % Match with ODE right-hand-side
        opti.subject_to(Pidot == h*f_x);
        opti.subject_to(0 == g_z);

        % State at end of collocation interval
        Xk_end = Z*D;
        %% New decision variable for state at end of a finite element
        Xk = opti.variable(n_x);
        Xs{end+1} = Xk;
        Xs_extended{end+1} = Xk;
        opti.subject_to(lbx <= Xk <= ubx);
        opti.set_initial(Xk, x0);
        ind_x= [ind_x,ind_total(end)+1:ind_total(end)+n_x];
        ind_total = [ind_total,ind_total(end)+1:ind_total(end)+n_x];
        % Continuity constraints
        opti.subject_to(Xk_end==Xk)
    end
end
%% Terminal Constraints
%  possible relaxation of terminal constraints

% terminal constraint for physical and numerical  time


%% Terminal Costs 
% basic terminal cost term

% Quadrature states

% Time-optimal problems



%% Elastic mode cost terms 

%% Regularization cost terms
% Regularization term for speed-of-time 
% idea: add sot to model reformulation and add f_q (or do it all here, to


% Cost term for grid regularization (step equilibration, heursitic step)

%% Collect all variables
Xs_extended = [Xs_extended{:}];
Xs = [Xs{:}];
Us = [Us{:}];
Zs = [Zs{:}];
%% CasADi functions for complementarity residuals (standard, cross_complementarity, joint)

%% Create NLP Solver instance 
opti.minimize(J);
opti.solver('ipopt');
sol = opti.solve();
%% Results
x_opt = sol.value(Xs);
x_opt_extended = sol.value(Xs_extended);
z_opt = sol.value(Zs);
u_opt = sol.value(Us);
J_opt = sol.value(J);

% Plot the solution
tgrid_x = linspace(0, T, N_FE*N_stg+1);
tgrid_u = linspace(0, T, N_stg+1);
clf;
hold on
plot(tgrid_x, x_opt(1,:), '--')
plot(tgrid_x, x_opt(2,:), '-')
stairs(tgrid_u, [u_opt nan], '-.')
xlabel('t')
legend('x1','x2','u')

%% Model: CasADi functions, dimenesions, auxilairy functions.

%% Solve initalization (bounds, inital guess, parameters)

%% Output of the function
varargout{1} = solver;
varargout{2} = solver_initalization;
varargout{3} = model;
varargout{4} = settings;
end