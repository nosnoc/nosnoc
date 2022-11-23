%
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
switching_case = 'sliding_mode'; 
%  Options: 'crossing' 'sliding_mode', 'spontaneous_switch' , 'leave_sliding_mode', 
%% NOSNOC settings
[settings] = default_settings_nosnoc();  %% Optionally call this function to have an overview of all options.
settings.n_s = 2;
settings.mpcc_mode = 3;
settings.homotopy_update_slope = 0.1;
settings.irk_scheme = 'Gauss-Legendre';
% settings.irk_scheme = 'Radau-IIA';
settings.irk_representation= 'differential';
settings.print_level = 2;
settings.step_equilibration = 1;
% discretization parameters
N_sim = 16;
T_sim = 1.5;

model.N_sim = N_sim;
model.N_finite_elements = 2;
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
        settings.use_previous_solution_as_initial_guess = 1;
        [results,stats] = integrator_fesd(model,settings);
        %
        figure
        plot(results.t_grid,results.x_res)
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
        [results,stats] = integrator_fesd(model,settings);
        %
        figure
        plot(results.t_grid,results.x_res)
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
        settings.irk_scheme = 'Explicit-RK';
        settings.n_s = 1;
        model.N_finite_elements = 3; % set 4, 5 for different outcomes
        settings.use_previous_solution_as_initial_guess = 1;
        [results,stats] = integrator_fesd(model,settings);
        %
        figure
        plot(results.t_grid,results.x_res)
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
        settings.use_previous_solution_as_initial_guess = 1;
        [results,stats] = integrator_fesd(model,settings);
        %
        figure
        plot(results.t_grid,results.x_res(1,:))
        grid on
        xlabel('$t$','Interpreter','latex')
        ylabel('$x(t)$','Interpreter','latex')
        grid on
    otherwise
        error('pick a value for example_num between 1 and 4.')
end



