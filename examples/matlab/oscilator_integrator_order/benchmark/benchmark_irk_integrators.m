%
%    This file is part of NOSNOC.
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
%
clc
%% Run experiments to test FESD integrator order for different IRK schemes
N_start = 11;
N_end = 1001;
N_samples = 5;



N_start = 11;
N_end = 115;
N_samples = 3;

 %% Gauss Legendre - Differential
% Gauss Legendre with FESD
 
legend_str = {'Midpoint Rule','Gauss-Legendre 4','Gauss-Legendre 6','Gauss-Legendre 8','Gauss-Legendre 10','Gauss-Legendre 12','Gauss-Legendre 14','Gauss-Legendre 16','Gauss-Legendre 18'};
settings.use_fesd = 1;
settings.print_level = 0;

settings.irk_scheme = 'Gauss-Legendre';
settings.irk_representation = 'differential';
settings.save_results = 1;
settings.n_s_vec = [1 2 3 4];
% settings.n_s_vec = [1];

settings.scenario_name = 'irk_gauss_legendre_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

%  Gauss Legendre without FESD
settings.use_fesd = 0;
settings.scenario_name = 'irk_gauss_legendre_std_differential';
[results] = integrator_order_experiment(settings,legend_str);
 

%% Radau II-A - Differential
% Radau II-A with FESD
legend_str = {'Implicit Euler','Radau-IIA 3','Radau-IIA 5','Radau-IIA 7','Radau-IIA 9','Radau-IIA 11','Radau-IIA 13','Radau-IIA 15','Radau-IIA 17'};
settings.use_fesd = 1;
settings.irk_scheme = 'Radau-IIA';
settings.irk_representation = 'differential';
settings.save_results = 1;
settings.n_s_vec = [1 2 3 4];
settings.scenario_name = 'irk_radau_iia_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
legend_str = {'Implicit Euler','Radau-IIA 3','Radau-IIA 5','Radau-IIA 7','Radau-IIA 9','Radau-IIA 11','Radau-IIA 13','Radau-IIA 15','Radau-IIA 17'};
[results] = integrator_order_experiment(settings,legend_str);

% Radau II-A without FESD
settings.use_fesd = 0;
settings.scenario_name = 'irk_radau_iia_std_differential';

[results] = integrator_order_experiment(settings,legend_str);
%  



%% Lobatto III-C (only differential)

 
legend_str = {'no Lobatto with n_s=1', 'Lobatto-IIIC 2','Lobatto-IIIC 4','Lobatto-IIIC 6','Lobatto-IIIC 8'};
settings.use_fesd = 1;
settings.irk_scheme = 'Lobatto-IIIC';
settings.irk_representation = 'differential';
settings.save_results = 1;
% settings.n_s_vec = [2 3 4 5];
settings.n_s_vec = [2 3 4 5];
settings.scenario_name = 'irk_Lobatto_iiic_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

%  Lobatto-IIIC without FESD
settings.use_fesd = 0;
settings.scenario_name = 'irk_Lobatto_iiic_std_differential';
% [results] = integrator_order_experiment(settings,legend_str);
 
%% Lobatto III-B (only differential)
 
legend_str = {'no Lobatto with n_s=1', 'Lobatto-IIIB 2','Lobatto-IIIB 4','Lobatto-IIIB 6','Lobatto-IIIB 8'};
settings.use_fesd = 1;
settings.irk_scheme = 'Lobatto-IIIB';
settings.irk_representation = 'differential';
settings.save_results = 1;
% settings.n_s_vec = [2 3 4 5];
settings.n_s_vec = [2 3 4 5];
settings.scenario_name = 'irk_Lobatto_iiib_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

%  Lobatto-IIIB without FESD
settings.use_fesd = 0;
settings.scenario_name = 'irk_Lobatto_iiib_std_differential';
[results] = integrator_order_experiment(settings,legend_str);
 


%% Lobatto-IIIA (only differential)
 
legend_str = {'no Lobatto with n_s=1', 'Lobatto-IIIA 2','Lobatto-IIIA 4','Lobatto-IIIA 6','Lobatto-IIIA 8'};
settings.use_fesd = 1;
settings.irk_scheme = 'Lobatto-IIIA';
settings.irk_representation = 'differential';
settings.save_results = 1;
settings.n_s_vec = [2 3 4 5];
settings.scenario_name = 'irk_Lobatto_iiia_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

%  Lobatto-IIIA without FESD
settings.use_fesd = 0;
settings.scenario_name = 'irk_Lobatto_iiia_std_differential';
[results] = integrator_order_experiment(settings,legend_str);
 


%% Lobatto-III (only differential)
 
legend_str = {'no Lobatto with n_s=1', 'Lobatto-III 2','Lobatto-III 4','Lobatto-III 6','Lobatto-III 8'};
settings.use_fesd = 1;
settings.irk_scheme = 'Lobatto-III';
settings.irk_representation = 'differential';
settings.save_results = 1;
settings.n_s_vec = [2 3 4 5];
% settings.n_s_vec = [2 3 4];
settings.scenario_name = 'irk_Lobatto_iii_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

%  Lobatto-III without FESD
settings.use_fesd = 0;
settings.scenario_name = 'irk_Lobatto_iii_std_differential';
[results] = integrator_order_experiment(settings,legend_str);
 


%% Radau-I (only differential)
% with fesd
 
legend_str = {'no Radau-I with n_s=1', 'Radau-I 3','Radau-I 5'};
settings.use_fesd = 1;
settings.irk_scheme = 'Radau-I';
settings.irk_representation = 'differential';
settings.save_results = 1;
settings.n_s_vec = [2 3];
settings.scenario_name = 'irk_radau_i_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

% without FESD
settings.use_fesd = 0;
settings.scenario_name = 'irk_radau_i_std_differential';
[results] = integrator_order_experiment(settings,legend_str);

%% Radau-IA
legend_str = {'Radau-IA 1', 'Radau-IA 3','Radau-IA 5'};
settings.use_fesd = 1;
settings.irk_scheme = 'Radau-IA';
settings.irk_representation = 'differential';
settings.save_results = 1;
settings.n_s_vec = [1 2 3];
settings.scenario_name = 'irk_radau_ia_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

% without FESD
settings.use_fesd = 0;
settings.scenario_name = 'irk_radau_ia_std_differential';
[results] = integrator_order_experiment(settings,legend_str);
 

%% Explicit RK (Experimental and for fun)
% with fesd
%  
legend_str = {'Explicit Euler', 'Heun 2','Kutta 3','Runge-Kutta 4','no','Nystrom 5'};
settings.use_fesd = 1;
settings.irk_scheme = 'Explicit-RK';
settings.irk_representation = 'differential';
settings.save_results = 1;
settings.n_s_vec = [1 2 3 4];
settings.scenario_name = 'erk_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

% without FESD
settings.use_fesd = 0;
settings.scenario_name = 'erk_std_differential';
[results] = integrator_order_experiment(settings,legend_str);
%  

% end


%% Radau II-A - Integral
% Radau II-A with FESD iintegral
%  
legend_str = {'Implicit Euler','Radau-IIA 3','Radau-IIA 5','Radau-IIA 7','Radau-IIA 9','Radau-IIA 11','Radau-IIA 13','Radau-IIA 15','Radau-IIA 17'};
settings.use_fesd = 1;
settings.irk_scheme = 'Radau-IIA';
settings.irk_representation = 'integral';
settings.save_results = 1;
settings.n_s_vec = [1 2 3 4];
settings.scenario_name = 'irk_radau_iia_fesd_integral';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
legend_str = {'Implicit Euler','Radau-IIA 3','Radau-IIA 5','Radau-IIA 7','Radau-IIA 9','Radau-IIA 11','Radau-IIA 13','Radau-IIA 15','Radau-IIA 17'};
[results] = integrator_order_experiment(settings,legend_str);

% Radau II-A without FESD integral
settings.use_fesd = 0;
settings.scenario_name = 'irk_radau_iia_std_integral';

% [results] = integrator_order_experiment(settings,legend_str);
%  
% end
%% Gauss Legendre - Integral
% Gauss Legendre with FESD
 
legend_str = {'Midpoint Rule','Gauss-Legendre 4','Gauss-Legendre 6','Gauss-Legendre 8','Gauss-Legendre 10','Gauss-Legendre 12','Gauss-Legendre 14','Gauss-Legendre 16','Gauss-Legendre 18'};
settings.use_fesd = 1;
settings.irk_scheme = 'Gauss-Legendre';
settings.irk_representation = 'integral';
settings.save_results = 1;
settings.n_s_vec = [1 2 3 4];
settings.scenario_name = 'irk_gauss_legendre_fesd_integral';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

%  Gauss Legendre without FESD
settings.use_fesd = 0;
settings.scenario_name = 'irk_gauss_legendre_std_integral';
[results] = integrator_order_experiment(settings,legend_str);

 
