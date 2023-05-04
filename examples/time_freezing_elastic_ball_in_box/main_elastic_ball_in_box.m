clc; close all; clear all;
import casadi.*
%% settings
% collocation settings
settings = NosnocOptions();
settings.n_s = 2;
settings.time_freezing_time_rescaling = 1;
settings.use_speed_of_time_variables =  1; 
settings.local_speed_of_time_variable = 1;  
settings.stagewise_clock_constraint = 1;
settings.time_freezing = 1;
settings.N_homotopy = 6;
% settings.homotopy_update_rule = 'superlinear';
%% Generate Model
% angulary velocity of reference
omega = -2*pi; % no impacts
% omega = -3*pi; % impacts
model = elastic_ball_in_box_model(omega);
%% Solve OCP via nosnoc
solver = NosnocSolver(model, settings);
[results,stats] = solver.solve();
%% Read and plot Result
unfold_struct(results,'base');
unfold_struct(model,'base');
plot_results_elastic_ball_in_box
