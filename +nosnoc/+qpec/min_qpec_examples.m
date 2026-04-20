clear; clc;
import casadi.*;

%% --- vertical QPEC dimensions ---
n0 = 50;     % free block
n1 = 50;      % first complementarity side
n2 = 50;      % second complementarity side
n_var  = n0 + n1 + n2;
n_ineq = 25;
n_comp = n1;   % complementarity pairs

%% --- PSD Q and random cost ---
R = randn(n_var);
Q = R'*R;
q = randn(n_var,1);

%% --- random linear constraints: lba <= A*x+b <= uba ---
A   = randn(n_ineq,n_var);
b   = randn(n_ineq,1);
lba = b - 1*abs(rand(n_ineq,1));
uba = b + 0*abs(rand(n_ineq,1));

%% --- vertical complementarity: g = x1, h = x2 ---
G = [zeros(n1,n0), eye(n1), zeros(n1,n2)];
H = [zeros(n2,n0+n1), eye(n2)];
g = zeros(n1,1);
h = zeros(n2,1);

%% ========== build base QPEC object ==========
qpec = nosnoc.qpec.Qpec(Q,A,G,H);
qpec.q   = q;
qpec.b   = b;
qpec.lba = lba;
qpec.uba = uba;
qpec.lbw = -10*ones(n_var,1);
qpec.ubw =  10*ones(n_var,1);
qpec.g   = g;
qpec.h   = h;
qpec.w0  = zeros(n_var,1);
qpec.y0  = zeros(n_comp,1);

methods = {'miqp','sos1','reg'};
RESULTS = struct();

for k = 1:numel(methods)
    mth = methods{k};

    opts = nosnoc.qpec.GurobiOptions();
    opts.method = mth;
    opts.bigM   = 10;
    opts.warmstart = true;
    opts.Nsqp = 6;
    opts.kappa = 0.01;
    opts.N_homotopy = 5;
    opts.lower_bounds_comps = true;
    opts.gurobi_params = struct('OutputFlag', 0);

    qpec.create_qpec_solver(opts,'gurobi');

    [res, st] = qpec.solve();

    x = res.x;
    y = res.y;
    obj = x'*Q*x + q'*x;
    comp_res = (G*x).*(H*x);

    RESULTS.(mth).x = x;
    RESULTS.(mth).y = y;
    RESULTS.(mth).obj = obj;
    RESULTS.(mth).comp_res = comp_res;
    RESULTS.(mth).cpu_time = st.cpu_time;    % added
end

%% ======== Print per-method details ========
fprintf('\n----- COMPARISON -----\n');
for k = 1:numel(methods)
    mth = methods{k};
    Rk = RESULTS.(mth);
    fprintf('\n### Method: %s ###\n', mth);
    fprintf('x*     = \n'); disp(Rk.x(:).');
    fprintf('y*     = \n'); disp(Rk.y(:).');
    fprintf('obj    = %g\n', Rk.obj);
    fprintf('comp   = (Gx).*(Hx) = \n'); disp(Rk.comp_res(:).');
    fprintf('cpu_time = %g seconds\n', Rk.cpu_time);
end

%% ======== Compact table comparison ========
mnames = methods;
nM = numel(mnames);

objs      = zeros(1,nM);
comp_norm = zeros(1,nM);
cpu_times = zeros(1,nM);

for j = 1:nM
    Rj = RESULTS.(mnames{j});
    objs(j)      = Rj.obj;
    comp_norm(j) = norm(Rj.comp_res, inf);
    cpu_times(j) = Rj.cpu_time;
end

fprintf('\n================ QPEC COMPARISON ================\n');
fprintf('%-15s', 'metric');
for j = 1:nM, fprintf('%15s', mnames{j}); end
fprintf('\n');

fprintf('%-15s', 'obj');
for j = 1:nM, fprintf('%15.6g', objs(j)); end
fprintf('\n');

fprintf('%-15s', '||comp||_inf');
for j = 1:nM, fprintf('%15.6g', comp_norm(j)); end
fprintf('\n');

fprintf('%-15s', 'cpu_time(s)');
for j = 1:nM, fprintf('%15.6g', cpu_times(j)); end
fprintf('\n');

fprintf('=================================================\n');

%% optional difference check
% norm(RESULTS.reg.x - RESULTS.miqp.x)
