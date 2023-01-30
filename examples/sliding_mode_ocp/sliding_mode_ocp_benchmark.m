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
%%
model = [];
[settings] = default_settings_nosnoc();  %% Optionally call this function to have an overview of all options.
settings.print_level = 3;
settings.mpcc_mode = 3;
settings.N_homotopy = 15;
settings.comp_tol = 1e-12;
settings.irk_representation = 'differential';
N_stages = 6;
settings.irk_scheme = 'Radau-IIA';
settings.n_s = 2;
scenario.N_finite_elements = 2;

settings.use_fesd = 1;

scenario.save_results = 1;
scenario.estimate_terminal_error = 1;


scenario.illustrate_regions  = 1;
scenario.terminal_constraint = 1;
scenario.linear_control = 1;
scenario.rho_v = 1;
scenario.rho_u = 1;
scenario.plot_results_sliding_mode = 0;

settings.step_equilibration  = 1;
settings.heuristic_step_equilibration = 0;


%% Basic
mpcc_mode_fesd = [5 5 5 3 3];
mpcc_mode_std = [3 3 3 3 3];
irk_schemes = {'Radau-IIA','Lobatto-IIIA','Lobatto-IIIC','Gauss-Legendre','Explicit-RK'};
n_s_vec = {[1:4],[2:4],[2:4],[1:4],[1 4]};
orders_vec = {n_s_vec{1}*2-1,n_s_vec{2}*2-2,n_s_vec{2}*2-2,n_s_vec{3}*2,n_s_vec{4}};
N_finite_elements_vec = 2:2:12;
N_experiments = 1:5;

%%
mpcc_mode_fesd = [5 5 3 3];
mpcc_mode_std = [3 3 3 3];
irk_schemes = {'Radau-IIA','Lobatto-IIIC','Gauss-Legendre','Explicit-RK'};
n_s_vec = {[1:4],[2:4],[1:4],[1 4]};
orders_vec = {n_s_vec{1}*2-1,n_s_vec{2}*2-2,n_s_vec{3}*2,n_s_vec{4}};
% N_finite_elements_vec = 2:2:12;
N_finite_elements_vec = 2:7;
N_experiments = 1:4;
%% for storing results
f_obj_fesd  = {};
cpu_fesd  = {};
comp_resiudal_fesd = {};
error_temrinal_fesd = {};
%
f_obj_std  = {};
cpu_std  = {};
comp_resiudal_std = {};
error_temrinal_std = {};

%% FESD
for k = 1:length(N_experiments)
    settings.irk_scheme = irk_schemes{k};

    f_obj_fesd_k  = [];
    cpu_fesd_k  = [];
    comp_resiudal_fesd_k = [];
    error_temrinal_fesd_k = [];
    % std
    f_obj_std_k  = [];
    cpu_std_k  = [];
    comp_resiudal_std_k = [];
    error_temrinal_std_k = [];
    settings.use_fesd = 1;
    settings.mpcc_mode= mpcc_mode_fesd(k);
    for jj = n_s_vec{k}
        f_obj_temp = [];
        cpu_temp = [];
        comp_resiudal_temp = [];
        error_temrinal_temp = [];
        for ii = N_finite_elements_vec
            settings.n_s = jj;
            scenario.N_finite_elements = ii;
            scenario.scenario_name = [settings.irk_scheme '_n_s_' num2str(settings.n_s) '_N_FE_' num2str(scenario.N_finite_elements) '_FESD_' num2str(settings.use_fesd)];

            [results] = sliding_mode_ocp_experiment(scenario,model,settings);
            f_obj_temp  = [f_obj_temp,results.f_opt];
            cpu_temp  = [cpu_temp,results.cpu_time];
            comp_resiudal_temp = [comp_resiudal_temp, results.comp_residual];
            error_temrinal_temp = [error_temrinal_temp, results.error_terminal];
        end
        f_obj_fesd_k  = [f_obj_fesd_k;f_obj_temp];
        cpu_fesd_k  = [cpu_fesd_k ;cpu_temp];
        comp_resiudal_fesd_k = [comp_resiudal_fesd_k ;comp_resiudal_temp];
        error_temrinal_fesd_k = [error_temrinal_fesd_k;error_temrinal_temp];
    end

    % std experiments
    settings.use_fesd = 0;
    settings.mpcc_mode= mpcc_mode_std(k);
    for jj = n_s_vec{k}
        f_obj_temp = [];
        cpu_temp = [];
        comp_resiudal_temp = [];
        error_temrinal_temp = [];
        for ii = N_finite_elements_vec
            settings.n_s = jj;
            scenario.N_finite_elements = ii;
            scenario.scenario_name = [settings.irk_scheme '_n_s_' num2str(settings.n_s) '_N_FE_' num2str(scenario.N_finite_elements) '_FESD_' num2str(settings.use_fesd)];
            [results] = sliding_mode_ocp_experiment(scenario,model,settings);
            f_obj_temp  = [f_obj_temp,results.f_opt];
            cpu_temp  = [cpu_temp,results.cpu_time];
            comp_resiudal_temp = [comp_resiudal_temp,results.comp_residual];
            error_temrinal_temp = [error_temrinal_temp, results.error_terminal];
        end
        f_obj_std_k  = [f_obj_std_k;f_obj_temp];
        cpu_std_k  = [cpu_std_k ;cpu_temp];
        comp_resiudal_std_k = [comp_resiudal_std_k ;comp_resiudal_temp];
        error_temrinal_std_k = [error_temrinal_std_k;error_temrinal_temp];
    end
    % store all data
    f_obj_fesd  = [f_obj_fesd;f_obj_fesd_k  ];
    cpu_fesd  = [cpu_fesd ;cpu_fesd_k  ];
    comp_resiudal_fesd = [comp_resiudal_fesd;comp_resiudal_fesd_k];
    error_temrinal_fesd = [error_temrinal_fesd;error_temrinal_fesd_k];

    f_obj_std  = [f_obj_std;f_obj_std_k];
    cpu_std  = [cpu_std;cpu_std_k];
    comp_resiudal_std = [comp_resiudal_std;comp_resiudal_std_k];
    error_temrinal_std = [error_temrinal_std;error_temrinal_std_k];
end
%% Plot results
sliding_moode_ocp_eval_results

%% Save results
% results std
results_benchmark.f_obj_std  = f_obj_std ;
results_benchmark.cpu_std  = cpu_std ;
results_benchmark.list_of_markers = list_of_markers;
results_benchmark.comp_resiudal_std  = comp_resiudal_std;
results_benchmark.error_temrinal_std  = error_temrinal_std;
% results fesd
results_benchmark.f_obj_fesd  = f_obj_fesd;
results_benchmark.cpu_fesd  = cpu_fesd ;
results_benchmark.comp_resiudal_fesd  = comp_resiudal_fesd;
results_benchmark.error_temrinal_fesd  = error_temrinal_fesd;
% scnearios
results_benchmark.irk_schemes = irk_schemes;
results_benchmark.N_experiments = N_experiments;
results_benchmark.N_finite_elements_vec  = N_finite_elements_vec ;
results_benchmark.n_s_vec  = n_s_vec ;

if scenario.save_results
    save(['results/results_benchmark_full_3.mat'],'results_benchmark')
    %     save([scenario_name '.mat'],'results')
end