function [results] = sliding_mode_ocp_experiment(scenario,model,settings)
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

import casadi.*

unfold_struct(scenario,'caller');
%% model equations
x_target = [-pi/6;-pi/3];
if linear_control
    v0  = [0;0];
else
    v0 = [];
end
model.x0 = [2*pi/3;pi/3;v0];
model.T = 4;
settings.N_stages = 6;
settings.N_finite_elements = N_finite_elements;
% Variable defintion
x1 = MX.sym('x1');
x2 = MX.sym('x2');

if linear_control
    v1 = MX.sym('v1');
    v2 = MX.sym('v2');
    x = [x1;x2;v1;v2];
else
    x = [x1;x2];
end
model.x = x;
% Control
u1 = MX.sym('u1');
u2 = MX.sym('u2');
model.u = [u1;u2];
if linear_control
    u_max = 5;
    model.lbu  = -u_max*ones(2,1);
    model.ubu  = u_max*ones(2,1);
else
    u_max = 2;
    model.lbu  = -u_max*ones(2,1);
    model.ubu  = u_max*ones(2,1);
end

% Switching Functions
% p = 2; a = -0.1; a1 = 0;
% b = -0.05; q = 3;

p = 2; a = 0.15; a1 = 0;
b = -0.05; q = 3;


c1 = x1+a*(x2-a1)^p;
c2 = x2+b*x1^q;

S1 = [1;-1];
S2 = [1;-1];
model.c = {c1,c2};
model.S = {S1,S2};

%% Modes of the ODEs layers (for all  i = 1,...,n_sys);
% part independet of the nonsmoothness
if linear_control
    f_11 = [-1+v1;0;u1;u2];
    f_12 = [1+v1;0;u1;u2];
    f_21 = [0;-1+v2;u1;u2];
    f_22 = [0;1+v2;u1;u2];
else
    f_11 = [-1+u1;0];
    f_12 = [1+u1;0];
    f_21 = [0;-1+u2];
    f_22 = [0;1+u2];
end
F1 = [f_11 f_12];
F2 = [f_21 f_22];
model.F = {F1,F2};

%% Objective
% model.f_q = u1^2+u2^2;
if linear_control
    model.f_q = rho_v*(v1^2+v2^2)+rho_u*(u1^2+u2^2);
else
    model.f_q = u1^2+u2^2;
end
if terminal_constraint
    model.g_terminal = [x(1:2)-x_target];
else
    model.f_q_T = 100*(x(1:2)-x_target(1:2))'*(x(1:2)-x_target(1:2));
end

% empty results
u_opt = [];
t_grid_optimizer = [];
x_res_optimizer = [];
error_terminal = [];
x_res_integrator = [];
t_grid_integrator = [];

%% Solve and plot
x0 = model.x0;
solver = NosnocSolver(model, settings);
[results,stats] = solver.solve();
f_opt = results.f_opt;
cpu_time = stats.cpu_time_total;

if plot_results_sliding_mode
    sliding_mode_ocp_plot
end
%% Establish integration error
% load results from optimal control problem
use_matlab_integrator = 1;
if estimate_terminal_error
    u_opt = results.u_opt;
    t_grid_optimizer = [results.t_grid];
    x_res_optimizer = [results.x_opt];
%     x_target = [-pi/6;-1];
    x_target = [-pi/6;-pi/3];

    if ~use_matlab_integrator
        model.T_sim = 4/6;
        model.N_sim = 32;
        settings.N_stages = 1;
        settings.N_finite_elements = 2;
        model.g_terminal = [];
        model.g_terminal_lb = [];
        model.g_terminal_ub = [];

        settings_integrator = settings;
        settings_integrator.mpcc_mode = 5;
        settings_integrator.comp_tol = 1e-15;
        settings_integrator.irk_scheme = IRKSchemes.RADAU_IIA;
        settings_integrator.n_s = 5;
        settings_integrator.use_fesd = 1;
        settings_integrator.print_level = 0;
        x_res_integrator = [];
        t_grid_integrator = [];
        t_end = 0;
        try
            for ii =  1:6
                model.lbu = u_opt(:,ii);
                model.ubu = u_opt(:,ii);
                model.u0 = u_opt(:,ii);
                [results_integrator,stats,model] = integrator_fesd(model,settings_integrator);
                model.x0 = results_integrator.x_res(:,end);
                x_res_integrator = [x_res_integrator,results_integrator.x_res];
                t_grid_integrator = [t_grid_integrator, results_integrator.t_grid+t_end];
                t_end = t_grid_integrator(end);
            end
            x_end = x_res_integrator(1:2,end);
            error_terminal = norm(x_target - x_end);
        catch
            try
                settings_integrator.mpcc_mode = 3;
                model.x0 = x0;
                for ii =  1:6
                    model.lbu = u_opt(:,ii);
                    model.ubu = u_opt(:,ii);
                    model.u0 = u_opt(:,ii);
                    [results_integrator,stats,model] = integrator_fesd(model,settings_integrator);
                    model.x0 = results_integrator.x_res(:,end);
                    x_res_integrator = [x_res_integrator,results_integrator.x_res];
                    t_grid_integrator = [t_grid_integrator, results_integrator.t_grid+t_end];
                    t_end = t_grid_integrator(end);
                end
                x_end = x_res_integrator(1:2,end);
                error_terminal = norm(x_target - x_end);
            catch
                error_terminal = nan;
            end
        end
    else

        tspan = [0 4/6];
        y0 = [2*pi/3;pi/3;0;0];
        sigma_int = 1e-11;
        tol = sigma_int/10;
        options = odeset('RelTol', tol, 'AbsTol', tol/10);
        x_res_integrator = [];
        t_grid_integrator = [];
        for ii = 1:6
            [t,y_res] = ode15s(@(t,y) ...
                ((0.5*(1+tanh((y(1)+a*(y(2)-a1)^p)/sigma_int))))*([-1+y(3);0;u_opt(1,ii);u_opt(2,ii)])+(1-(0.5*(1+tanh((y(1)+a*(y(2)-a1)^p)/sigma_int))))*([1+y(3);0;u_opt(1,ii);u_opt(2,ii)])+...
                ((0.5*(1+tanh((y(2)+b*y(1)^q)/sigma_int))))*([0;-1+y(4);u_opt(1,ii);u_opt(2,ii)])+(1-(0.5*(1+tanh((y(2)+b*y(1)^q)/sigma_int))))*([0;1+y(4);u_opt(1,ii);u_opt(2,ii)])...
                ,tspan, y0,options);

            y0 = y_res(end,:)';
            x_res_integrator = [x_res_integrator,y_res'];
            t_grid_integrator = [t_grid_integrator,t'+(ii-1)*4/6];
        end
        x_end = x_res_integrator(1:2,end);
        error_terminal = norm(x_target - x_end);
    end
    fprintf('Error %2.4e \n',error_terminal);

end

%%
fprintf('Objective value %2.4f \n',f_opt);

%% Save results (cpu time, objective value, comp tolerance)
results.cpu_time = cpu_time;
results.comp_residual = stats.complementarity_stats(end);
results.u_opt = u_opt;
results.f_opt = f_opt;
results.t_grid_optimizer = t_grid_optimizer;
results.x_res_optimizer = x_res_optimizer;
results.error_terminal  = error_terminal;
results.x_res_integrator = x_res_integrator;
results.t_grid_integrator = t_grid_integrator;

if save_results
    save(['results/' scenario_name '.mat'],'results')
    %     save([scenario_name '.mat'],'results')
end
end

