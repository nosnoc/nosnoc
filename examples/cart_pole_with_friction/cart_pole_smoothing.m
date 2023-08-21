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

%% model stuff
% % Fixed friction force
% F_friction = 2;
% c = v(1);
% xdot = [v;inv(M)*(f_all - [F_friction*tanh(c/sigma);0])];
% 

sigma = SX.sym('sigma',1);

model = get_cart_pole_with_friction_model();
coeff_0 = (1+tanh(model.c/sigma)) / 2;
coeff_1 = 1 - coeff_0;

model = struct(model);
model.xdot_smoothened = coeff_0 * model.F(:, 1) + coeff_1 * model.F(:, 2);

%%
n_s = 2;
N = 30; % number of control intervals
x_ref = [0; 180/180*pi; 0; 0]; % target position

%% NLP formulation & solver creation
nlp = setup_collocation_nlp(model, sigma, n_s, N);

casadi_nlp = struct('f', nlp.f, 'x', nlp.w, 'p', nlp.p, 'g', nlp.g);

solver = nlpsol('solver', 'ipopt', casadi_nlp);

%% solve nlp
sigma0 = 1e-3;
% Solve the NLP
sol = solver('x0', nlp.w0, 'lbx', nlp.lbw, 'ubx', nlp.ubw,...
    'lbg', nlp.lbg, 'ubg', nlp.ubg, 'p', sigma0);
w_opt = full(sol.x);

%% extract results
n_x = length(model.x);
n_u = length(model.u);
q1_opt = w_opt(1:(n_x+n_u)+n_x*(n_s):end)';
q2_opt = w_opt(2:(n_x+n_u)+n_x*(n_s):end)';
v1_opt = w_opt(3:(n_x+n_u)+n_x*(n_s):end)';
v2_opt = w_opt(4:(n_x+n_u)+n_x*(n_s):end)';
u_opt = w_opt(5:(n_x+n_u)+n_x*(n_s):end)';

results = struct();
results.x = [q1_opt; q2_opt; v1_opt; v2_opt];
results.t_grid = linspace(0, model.T, N+1);
results.t_grid_u = linspace(0, model.T, N+1);
results.u = u_opt;

model.h_k = model.T / (N * n_s);
plot_cart_pole_trajecetory(results, model, x_ref);

%%
% figure
% F_friction_opt = F_friction*tanh(v1_opt/sigma0);
% plot(t_grid,F_friction_opt)
% hold on
% stairs(t_grid,[u_opt, nan]');
% 
% grid on
% yline(0,'k--');
% ylabel('$F_{\mathrm{fr}}(t)$','Interpreter','latex')
% xlabel('$t$','Interpreter','latex')

%% objective plot
% f_opt = [];
% sigma_vec = [];
% sigma0 = 1e1;
% kappa = 0.25;
% w0_def = w0;
% for ii = 1:15
%     sigma_vec = [sigma_vec, sigma0];
%     sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg, 'p', sigma0);
%     f_opt = [f_opt, full(sol.f)];
%     sigma0 = kappa*sigma0;
%     % w0 = full(sol.x);
% end
% figure
% semilogx(sigma_vec,f_opt)
% grid on



function nlp = setup_collocation_nlp(model, sigma, n_s, N)
    import casadi.*
    %% collocation using casadi
    % Get collocation points
    tau_root = [0 collocation_points(n_s, 'radau')];
    
    % Coefficients of the collocation equation
    C = zeros(n_s+1,n_s+1);
    
    % Coefficients of the continuity equation
    D = zeros(n_s+1, 1);
    
    % Coefficients of the quadrature function
    B = zeros(n_s+1, 1);
    
    % Construct polynomial basis
    for j=1:n_s+1
        % Construct Lagrange polynomials to get the polynomial basis at the collocation point
        coeff = 1;
        for r=1:n_s+1
            if r ~= j
                coeff = conv(coeff, [1, -tau_root(r)]);
                coeff = coeff / (tau_root(j)-tau_root(r));
            end
        end
        % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
        D(j) = polyval(coeff, 1.0);
    
        % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
        pder = polyder(coeff);
        for r=1:n_s+1
            C(j,r) = polyval(pder, tau_root(r));
        end
    
        % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
        pint = polyint(coeff);
        B(j) = polyval(pint, 1.0);
    end
    
    
    %% Extract from model
    L = model.f_q;
    n_x = length(model.x);
    n_u = length(model.u);

    % Continuous time dynamics
    f = Function('f', {model.x, model.u, sigma}, {model.xdot_smoothened, L});
    f_terminal = Function('f', {model.x,}, {model.f_q_T});

    % Control discretization
    h = model.T/N;
    x0 = model.x0;
    
    %% Casadi NLP formulation
    w = {};
    w0 = [];
    lbw = [];
    ubw = [];
    J = 0;
    g = {};
    lbg = [];
    ubg = [];

    % "Lift" initial conditions
    Xk = SX.sym('X0', n_x);
    w = {w{:}, Xk};
    lbw = [lbw; x0];
    ubw = [ubw; x0];
    w0 = [w0; x0];
    
    % Formulate the NLP
    for k=0:N-1
        % New NLP variable for the control
        Uk = SX.sym(['U_' num2str(k)],n_u);
        w = {w{:}, Uk};
        lbw = [lbw; model.lbu];
        ubw = [ubw; model.ubu];
        w0 = [w0;  0];
    
        % State at collocation points
        Xkj = {};
        for j=1:n_s
            Xkj{j} = SX.sym(['X_' num2str(k) '_' num2str(j)], n_x);
            w = {w{:}, Xkj{j}};
            lbw = [lbw; model.lbx];
            ubw = [ubw; model.ubx];
            w0 = [w0; x0];
        end
    
        % Loop over collocation points
        Xk_end = D(1)*Xk;
        for j=1:n_s
            % Expression for the state derivative at the collocation point
            xp = C(1,j+1)*Xk;
            for r=1:n_s
                xp = xp + C(r+1,j+1)*Xkj{r};
            end
    
            % Append collocation equations
            [fj, qj] = f(Xkj{j},Uk, sigma);
            g = {g{:}, h*fj - xp};
            lbg = [lbg; zeros(n_x,1)];
            ubg = [ubg; zeros(n_x,1)];
    
            % Add contribution to the end state
            Xk_end = Xk_end + D(j+1)*Xkj{j};
    
            % Add contribution to quadrature function
            J = J + B(j+1)*qj*h;
        end
    
        % New NLP variable for state at end of interval
        Xk = SX.sym(['X_' num2str(k+1)], n_x);
        w = {w{:}, Xk};
        lbw = [lbw; model.lbx];
        ubw = [ubw; model.ubx];
        w0 = [w0; x0];
    
        % Add equality constraint
        g = {g{:}, Xk_end-Xk};
        lbg = [lbg; zeros(n_x,1)];
        ubg = [ubg; zeros(n_x,1)];
    end

    J = J + f_terminal(Xk_end);
    
    nlp = struct();
    nlp.f = J;
    nlp.w = vertcat(w{:});
    nlp.w0 = w0;
    nlp.lbw = lbw;
    nlp.ubw = ubw;
    nlp.p = sigma;
    nlp.g = vertcat(g{:});
    nlp.lbg = lbg;
    nlp.ubg = ubg;
end
