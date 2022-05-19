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
function [results] = local_minima_experiment(scenario,settings,model)
import casadi.*
unfold_struct(scenario,'caller');

%% Local min expriment
% solve the ocp for different inital guesses and see does it get stuck due
% to wrong derivatives
%% Convergence test:
L_numeric = [];
complementarity_stats = [];

[solver,solver_initalization, model,settings] = create_nlp_nosnoc(model,settings);
unfold_struct(model,'caller');
unfold_struct(settings,'caller');
unfold_struct(solver_initalization,'caller');

x0_star = [];
solver_initalization.lbw(1) =  -inf;
solver_initalization.ubw(1) =  inf;

for jj = 1:N_samples
    x0 = x0_vec(jj);
    solver_initalization.w0(ind_x) = x0;
    solver_initalization.lbw(1) =  x0;
    solver_initalization.ubw(1) =  x0;
    % forward simulation for initalization
    [sol,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);

    % solve NLP
    solver_initalization.w0 = full(sol.x);
    solver_initalization.lbw(1) =  -inf;
    solver_initalization.ubw(1) =  inf;
    [sol,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);
    w_opt = full(sol.x);
    x1_opt = w_opt(ind_x);

    x0_star = [x0_star;x1_opt(1)];
    try
    f_opt = full(J_fun(w_opt));
    catch
        f_opt  = 0;
    end
    L_numeric = [L_numeric;f_opt];
    complementarity_iter = full(comp_res(w_opt));
    complementarity_stats = [complementarity_stats;complementarity_iter];
    ts = -x0/3;
    L_analytic =(8/3*ts^3 +8/3*x0*ts^2 + 8/9*x0^2*ts +1/3*T^3+1/3*x0*T^2+x0^2*T/9) + (T+(x0-5)/3)^2;
end
x0_analytic = -1.427572232082767;
figure
plot(x0_vec,x0_analytic*ones(1,N_samples),'LineWidth',1.0);
hold on
plot(x0_vec,x0_star,'LineWidth',1.0);
axis equal
xlabel('$x_0$','interpreter','latex');
ylabel('$x_0^*$','interpreter','latex');
if settings.use_fesd
    legend({'Analytic Solution','FESD Solution'},'interpreter','latex');
else
    legend({'Analytic Solution','Standard Solution'},'interpreter','latex');
end
grid on

%%
results.x0_vec = x0_vec;
results.x0_analytic = x0_analytic*ones(1,N_samples);
results.x0_star = x0_star;
results.L_numeric= L_numeric;
results.L_analytic= L_analytic*ones(1,N_samples);
results.complementarity_stats= complementarity_stats;


if save_results
    %     save(['results/' scenario_name '.mat'],'results')
    save([scenario_name '.mat'],'results')
end

end

