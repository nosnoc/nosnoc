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
%
%% Run experiments to test FESD integrator order for different IRK schemes
N_start = 11;
% N_end = 2001;
N_end = 1021;
N_samples = 5;

%% Gauss Legendre - Differential
% Gauss Legendre with FESD
try
legend_str = {'Midpoint Rule','Gauss-Legendre 4','Gauss-Legendre 6','Gauss-Legendre 8','Gauss-Legendre 10','Gauss-Legendre 12','Gauss-Legendre 14','Gauss-Legendre 16','Gauss-Legendre 18'};
settings.use_fesd = 1;
settings.irk_scheme = 'Gauss-Legendre';
settings.irk_representation = 'differential';
settings.save_results = 1;
settings.n_s_vec = [1 2 3 4];
settings.scenario_name = 'irk_gauss_legendre_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

%  Gauss Legendre without FESD
settings.use_fesd = 0;
settings.scenario_name = 'irk_gauss_legendre_std_differential';
[results] = integrator_order_experiment(settings,legend_str);
catch
end

%% Radau II-A - Differential
% Radau II-A with FESD
try
legend_str = {'Implicit Euler','Radau-IIA 3','Radau-IIA 5','Radau-IIA 7','Radau-IIA 9','Radau-IIA 11','Radau-IIA 13','Radau-IIA 15','Radau-IIA 17'};
settings.use_fesd = 1;
settings.irk_scheme = 'Radau-IIA';
settings.irk_representation = 'differential';
settings.save_results = 1;
% settings.n_s_vec = [1 2 3 4 5 6];
settings.n_s_vec = [1 2 3 4 5];
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
catch
end

%% Radau-I (only differential)
% with fesd
try
legend_str = {'no Radau-I with n_s=1', 'Radau-I 3','Radau-I 5'};
settings.use_fesd = 1;
settings.irk_scheme = 'Radua-I';
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
catch

end


%% Radau-IA (only differential)
% with fesd
try
legend_str = {'Radau-IA 1', 'Radau-IA 3','Radau-IA 5'};
settings.use_fesd = 1;
settings.irk_scheme = 'Radau-IA';
settings.irk_representation = 'differential';
settings.save_results = 1;
settings.n_s_vec = [1 2 3];
settings.n_s_vec = [2];
settings.scenario_name = 'irk_radau_ia_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

% without FESD
settings.use_fesd = 0;
settings.scenario_name = 'irk_radau_ia_std_differential';
[results] = integrator_order_experiment(settings,legend_str);
catch

end


%% Lobbato III-C (only differential)

try
legend_str = {'no Lobbato with n_s=1', 'Lobbato-IIIC 2','Lobbato-IIIC 4','Lobbato-IIIC 6','Lobbato-IIIC 8'};
settings.use_fesd = 1;
settings.irk_scheme = 'Lobbato-IIIC';
settings.irk_representation = 'differential';
settings.save_results = 1;
% settings.n_s_vec = [2 3 4 5];
settings.n_s_vec = [2 3 4 5];
settings.scenario_name = 'irk_lobbato_iiic_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

%  Lobbato-IIIC without FESD
settings.use_fesd = 0;
settings.scenario_name = 'irk_lobbato_iiic_std_differential';
% [results] = integrator_order_experiment(settings,legend_str);
catch

end

%% Lobbato III-B (only differential)
try
legend_str = {'no Lobbato with n_s=1', 'Lobbato-IIIB 2','Lobbato-IIIB 4','Lobbato-IIIB 6','Lobbato-IIIB 8'};
settings.use_fesd = 1;
settings.irk_scheme = 'Lobbato-IIIB';
settings.irk_representation = 'differential';
settings.save_results = 1;
% settings.n_s_vec = [2 3 4 5];
settings.n_s_vec = [2 3 4 5];
settings.scenario_name = 'irk_lobbato_iiib_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

%  Lobbato-IIIB without FESD
settings.use_fesd = 0;
settings.scenario_name = 'irk_lobbato_iiib_std_differential';
[results] = integrator_order_experiment(settings,legend_str);
catch

end
%% Lobbato-IIIA (only differential)
try
legend_str = {'no Lobbato with n_s=1', 'Lobbato-IIIA 2','Lobbato-IIIA 4','Lobbato-IIIA 6','Lobbato-IIIA 8'};
settings.use_fesd = 1;
settings.irk_scheme = 'Lobbato-IIIA';
settings.irk_representation = 'differential';
settings.save_results = 1;
settings.n_s_vec = [2 3 4 5];
settings.scenario_name = 'irk_lobbato_iiia_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

%  Lobbato-IIIA without FESD
settings.use_fesd = 0;
settings.scenario_name = 'irk_lobbato_iiia_std_differential';
[results] = integrator_order_experiment(settings,legend_str);
catch

end
%% Lobbato-III (only differential)
try
legend_str = {'no Lobbato with n_s=1', 'Lobbato-III 2','Lobbato-III 4','Lobbato-III 6','Lobbato-III 8'};
settings.use_fesd = 1;
settings.irk_scheme = 'Lobbato-III';
settings.irk_representation = 'differential';
settings.save_results = 1;
settings.n_s_vec = [2 3 4 5];
% settings.n_s_vec = [2 3 4];
settings.scenario_name = 'irk_lobbato_iii_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

%  Lobbato-III without FESD
settings.use_fesd = 0;
settings.scenario_name = 'irk_lobbato_iii_std_differential';
[results] = integrator_order_experiment(settings,legend_str);
catch

end



%% Explicit RK (Experimental and for fun)
% with fesd
try
legend_str = {'Explicit Euler', 'Heun 2','Kutta 3','Runge-Kutta 4','Nystrom 5'};
settings.use_fesd = 1;
settings.irk_scheme = 'Explicit-RK';
settings.irk_representation = 'differential';
settings.save_results = 1;
settings.n_s_vec = [1 2 3 4 6];
settings.scenario_name = 'erk_fesd_differential';
settings.N_start = N_start;
settings.N_end = N_end;
settings.N_samples = N_samples;
[results] = integrator_order_experiment(settings,legend_str);

% without FESD
settings.use_fesd = 0;
settings.scenario_name = 'erk_std_differential';
[results] = integrator_order_experiment(settings,legend_str);
catch

end


%% Radau II-A - Integral
% Radau II-A with FESD iintegral
try
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
catch
end
%% Gauss Legendre - Integral
% Gauss Legendre with FESD
try
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

catch

end