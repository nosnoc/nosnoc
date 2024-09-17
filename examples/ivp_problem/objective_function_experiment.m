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

function [results] = objective_function_experiment(scenario, problem_options, solver_options, model)
import casadi.*
unfold_struct(scenario,'caller');
%% objective experimetn
L_numeric = [];
complementarity_stats = [];

error_state = [];
error_objective = [];

ocp_solver = nosnoc.ocp.Solver(model, problem_options, solver_options);

for jj = 1:N_samples
    x0 = x0_vec(jj);
    ocp_solver.set('x', 'init', {0,0,problem_options.n_s}, x0);
    ocp_solver.set('x', 'lb', {0,0,problem_options.n_s}, x0);
    ocp_solver.set('x', 'ub', {0,0,problem_options.n_s}, x0);
    ocp_solver.solve();
    
    x1_opt = ocp_solver.get('x');
    
    f_opt = ocp_solver.get_objective();
    L_numeric = [L_numeric;f_opt];
    complementarity_iter = ocp_solver.stats.complementarity_stats(end);
    complementarity_stats = [complementarity_stats;complementarity_iter];
    T = problem_options.T;
    ts = -x0/3;
    L_analytic =(8/3*ts^3 +8/3*x0*ts^2 + 8/9*x0^2*ts +1/3*T^3+1/3*x0*T^2+x0^2*T/9) + (T+(x0-5)/3)^2;
    error_state = [error_state;abs(x1_opt(end)-(T-ts))];
    error_objective = [error_objective;abs(f_opt-L_analytic)];
end

figure
subplot(211)
semilogy(x0_vec,error_state)
grid on
ylabel('Error state','interpreter','latex');
subplot(212)
semilogy(x0_vec,error_objective)
grid on
xlabel('$x_0$','interpreter','latex');
ylabel('Error objective','interpreter','latex');
%
figure
plot(x0_vec,L_numeric,'LineWidth',1.0);
ts = -x0_vec/3;
L_analytic =(8/3*ts.^3 +8/3*x0_vec.*ts.^2 + 8/9*x0_vec.^2.*ts +1/3*T^3+1/3*x0_vec.*T^2+x0_vec.^2*T/9) + (T+(x0_vec-5)/3).^2;
hold on
plot(x0_vec,L_analytic,'k--','LineWidth',2.5);
xlabel('$x_0$','interpreter','latex');
ylabel('Objective','interpreter','latex');
if problem_options.use_fesd
    legend({'$V_{\mathrm{FESD}}(x_0)$','$V^{*}(x_0)$'},'interpreter','latex');
else
    legend({'$V_{\mathrm{Num}}(x_0)$','$V^{*}(x_0)$'},'interpreter','latex');
end
grid on

%% save results
results.x0_vec = x0_vec;
results.L_analytic = L_analytic;
results.L_numeric = L_numeric;
results.error_state= error_state;
results.error_objective= error_objective;
results.error_state= error_state;
results.complementarity_stats= complementarity_stats;


if save_results
%     save(['results/' scenario_name '.mat'],'results')
    save([scenario_name '.mat'],'results')
end
end

