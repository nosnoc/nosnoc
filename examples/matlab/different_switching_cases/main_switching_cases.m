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
clear all
clc
close all
import casadi.*
%% Info
% The goal of this example is to ilustrate the 4 different cases of switching that can occur in a PSS, on some simple examples.
% by chossing example_num from 1 to 4 below one can see a case of
% 1) crossing a disconituity
% 2) sliding mode
% 3) sliding on a surfce of disconinuity where a spontenus switch can happen (nonuqnie solutions)
% 4) unique leaving of a sliding mode
example_num = 4;
%% NOS-NOC settings
[settings] = default_settings_fesd();  %% Optionally call this function to have an overview of all options.
settings.n_s = 3;
settings.mpcc_mode = 3;
settings.kappa = 0.1;
% settings.use_fesd = 0;
% settings.irk_scheme = 'Lobatto-IIIC';
settings.irk_scheme = 'Gauss-Legendre';
settings.irk_representation= 'differential';
settings.print_level = 2;
% discretization parameters
N_sim = 16;
N_stages = 2;
T_sim = 1.5;

% model.T_sim = T_sim ;
model.N_sim = N_sim;
model.N_stages = N_stages;
model.N_finite_elements = 1;
model.T_sim = T_sim;

switch example_num
    case 1
        %% Crossing a disconitnuity
        

        model.x0 = [-1];
        x = MX.sym('x',1);
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
    case 2
        %% Sliding mode
        
        model.x0 = [-0.5];
        x = MX.sym('x',1);
        model.x = x;
        model.c = x;
        model.S = [-1; 1];
        f_1 = [1]; f_2 = [-1];
        model.F = [f_1 f_2];

%         settings.use_previous_solution_as_initial_guess = 1;
        [results,stats] = integrator_fesd(model,settings);
        %
        figure
        plot(results.t_grid,results.x_res)
        grid on
        xlabel('$t$','Interpreter','latex')
        ylabel('$x(t)$','Interpreter','latex')
        grid on
    case 3
        %% spontenus switch
        
        model.x0 = [0.0];
        x = MX.sym('x',1);
        model.x = x;
        model.c = x;
        model.S = [-1; 1];
        f_1 = [-1]; f_2 = [1];
        model.F = [f_1 f_2];
        settings.use_previous_solution_as_initial_guess = 0;
        [results,stats] = integrator_fesd(model,settings);
        %
        figure
        plot(results.t_grid,results.x_res)
        grid on
        xlabel('$t$','Interpreter','latex')
        ylabel('$x(t)$','Interpreter','latex')
        ylim([-1 1])
        grid on

    case 4
        %% leaving sliding mode in a unique way
        model.x0 = [0;0];
        x = MX.sym('x',1);
        t = MX.sym('t',1);
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



