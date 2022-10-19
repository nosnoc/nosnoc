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
%
function [results] = objective_function_experiment(scenario,settings,model)
import casadi.*
unfold_struct(scenario,'caller');
%% objective experimetn
L_numeric = [];
complementarity_stats = [];

error_state = [];
error_objective = [];

[solver,solver_initialization, model,settings] = create_nlp_nosnoc(model,settings);
unfold_struct(model,'caller');
unfold_struct(settings,'caller');
unfold_struct(solver_initialization,'caller');

for jj = 1:N_samples
    x0 = x0_vec(jj);
    solver_initialization.lbw(1) =  x0;
    solver_initialization.ubw(1) =  x0;
    % solve NLP
    [sol,stats,solver_initialization] = homotopy_solver(solver,model,settings,solver_initialization);
    w_opt = full(sol.x);
    x1_opt = w_opt(ind_x);
    solver_initialization.w0 = w_opt;
    f_opt = full(J_fun(w_opt));
    L_numeric = [L_numeric;f_opt];
    complementarity_iter = full(comp_res(w_opt));
    complementarity_stats = [complementarity_stats;complementarity_iter ];
    ts = -x0/3;
    L_analytic =(8/3*ts^3 +8/3*x0*ts^2 + 8/9*x0^2*ts +1/3*T^3+1/3*x0*T^2+x0^2*T/9) + (T+(x0-5)/3)^2;
    error_state = [error_state;abs(x1_opt(end)-(T-ts)) ];
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
if settings.use_fesd
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

