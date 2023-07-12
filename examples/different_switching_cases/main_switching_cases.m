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

clear all
clc
close all
import casadi.*
%% Info
% The goal of this example is to illustrate the 4 different cases of switching that can occur in a PSS, on some simple examples.
% by chosing example_num from 1 to 4 below one can see a case of
% 1) crossing a disconituity
% 2) sliding mode
% 3) sliding on a surfce of disconinuity where a spontaneous switch can happen (nonuqnie solutions)
% 4) unique leaving of a sliding mode
switching_case = 'crossing'; 
%  Options: 'crossing' 'sliding_mode', 'spontaneous_switch' , 'leave_sliding_mode', 
%% NOSNOC settings
problem_options = NosnocProblemOptions();
solver_options = NosnocSolverOptions();
model = NosnocModel();

problem_options.n_s = 2;
solver_options.homotopy_update_slope = 0.1;
problem_options.irk_scheme = IRKSchemes.GAUSS_LEGENDRE;
% settings.irk_scheme = IRKSchemes.RADAU_IIA;
problem_options.irk_representation= 'differential';
solver_options.print_level = 3;
solver_options.store_integrator_step_results = 1;
% discretization parameters
N_sim = 16;
T_sim = 1.5;


model.N_sim = N_sim;
problem_options.N_finite_elements = 2;
model.T_sim = T_sim;

switch switching_case
    case 'crossing'
        %% Crossing a discontinuity
        model.x0 = -1;
        x = SX.sym('x',1);
        model.x = x;
        model.c = x;
        model.S = [-1; 1];
        f_1 = [2]; f_2 = [0.2];
        model.F = [f_1 f_2];
        solver_options.use_previous_solution_as_initial_guess = 1;
        integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
        [results,stats] = integrator.solve();
        %
        figure
        plot(results.t_grid,results.x)
        grid on
        xlabel('$t$','Interpreter','latex')
        ylabel('$x(t)$','Interpreter','latex')
        grid on
    case 'sliding_mode'
        %% Sliding mode
        model.x0 = -0.5;
        x = SX.sym('x',1);
        model.x = x;
        model.c = x;
        model.S = [-1; 1];
        f_1 = [1]; f_2 = [-1];
        model.F = [f_1 f_2];
        integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
        [results,stats] = integrator.solve();
        %
        figure
        plot(results.t_grid,results.x)
        grid on
        xlabel('$t$','Interpreter','latex')
        ylabel('$x(t)$','Interpreter','latex')
        grid on
    case 'spontaneous_switch'
        %% spontaneous switch
        model.x0 = 0;
        x = SX.sym('x',1);
        model.x = x;
        model.c = x;
        model.S = [-1; 1];
        f_1 = [-1]; f_2 = [1];
        model.F = [f_1 f_2];
        % implcit methods more accurate, explicit Euler enables "random"
        % leaving
        problem_options.irk_scheme = 'EXPLICIT_RK';
        problem_options.n_s = 1;
        problem_options.N_finite_elements = 3; % set 4, 5 for different outcomes
        solver_options.use_previous_solution_as_initial_guess = 1;
        integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
        [results,stats] = integrator.solve();
        %
        figure
        plot(results.t_grid,results.x)
        grid on
        xlabel('$t$','Interpreter','latex')
        ylabel('$x(t)$','Interpreter','latex')
        ylim([-1 1])
        grid on

    case 'leave_sliding_mode'
        %% leaving sliding mode in a unique way
        model.x0 = [0;0];
        x = SX.sym('x',1);
        t = SX.sym('t',1);
        model.x = [x;t];
        model.c = x;
        model.S = [-1; 1];
        f_1 = [1+t;1]; f_2 = [-1+t;1];
        model.F = [f_1 f_2];
        solver_options.use_previous_solution_as_initial_guess = 1;
        integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
        [results,stats] = integrator.solve();
        %
        figure
        plot(results.t_grid,results.x(1,:))
        grid on
        xlabel('$t$','Interpreter','latex')
        ylabel('$x(t)$','Interpreter','latex')
        grid on
    otherwise
        error('pick a value for example_num between 1 and 4.')
end



