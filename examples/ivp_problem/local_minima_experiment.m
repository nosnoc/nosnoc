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

[solver, solver_initialization, model, settings] = create_nlp_nosnoc(model,settings);
unfold_struct(model,'caller');
unfold_struct(settings,'caller');
unfold_struct(solver_initialization,'caller');

x0_star = [];
solver_initialization.lbw(1) =  -inf;
solver_initialization.ubw(1) =  inf;

for jj = 1:N_samples
    x0 = x0_vec(jj);
    solver_initialization.w0(ind_x) = x0;
    solver_initialization.lbw(1) =  x0;
    solver_initialization.ubw(1) =  x0;
    % forward simulation for initialization
    [sol,stats,solver_initialization] = homotopy_solver(solver,model,settings,solver_initialization);

    % solve NLP
    solver_initialization.w0 = full(sol.x);
    solver_initialization.lbw(1) =  -inf;
    solver_initialization.ubw(1) =  inf;
    [sol,stats,solver_initialization] = homotopy_solver(solver,model,settings,solver_initialization);
    w_opt = full(sol.x);
    x1_opt = w_opt(ind_x);

    x0_star = [x0_star;x1_opt(1)];
    try
    f_opt = full(J_fun(w_opt));
    catch
        f_opt  = 0;
    end
    L_numeric = [L_numeric;f_opt];
    complementarity_iter = full(comp_res(w_opt,[]));
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

