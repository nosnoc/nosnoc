function  [f_lin_opt, d_lpcc, y_lpcc, lpcc_solution_exists, lpcc_cpu_time] = solve_lpcc(x_k,...
            dims,...
            d_lpcc, y_lpcc, Delta_TR_k_l,...
            lbx, ubx,...
            nabla_f_k, A_k, a_k, B_k, b_k, B_comp_k, b_comp_k)


sense = [repmat('=',1,dims.n_eq),repmat('>',1,dims.n_ineq),repmat('>',1,2*dims.n_comp)];
vtype = [repmat('C',1,dims.n_primal),repmat('B',1,dims.n_comp)];
vtype_num = [repmat(0,1,dims.n_primal),repmat(1,1,dims.n_comp)];


% lb
lb = lbx;
lb(dims.ind_x0) = max([lbx(dims.ind_x0)-x_k(dims.ind_x0),-Delta_TR_k_l*ones(size(dims.ind_x0))]');
lb(dims.ind_x1) = max([lbx(dims.ind_x1)-x_k(dims.ind_x1),-Delta_TR_k_l*ones(dims.n_comp,1),-x_k(dims.ind_x1)],[],2);
lb(dims.ind_x2) = max([lbx(dims.ind_x2)-x_k(dims.ind_x2),-Delta_TR_k_l*ones(dims.n_comp,1),-x_k(dims.ind_x2)],[],2);
lb = [lb;-inf*ones(dims.n_comp,1)];
% ub
ub = [min(ubx-x_k,Delta_TR_k_l);inf*ones(dims.n_comp,1)];
%% set up gurobi problem
model.A = sparse([[A_k,zeros(dims.n_eq,dims.n_comp)];[B_k,zeros(dims.n_ineq,dims.n_comp)];B_comp_k]);
model.sense = sense;
model.rhs = [-a_k; -b_k; b_comp_k];
model.obj = [nabla_f_k;zeros(dims.n_comp,1)];
model.vtype = vtype;
model.modelsense = 'min';
model.lb = lb;
model.ub = ub;
params.outputflag = 0;
% params.BarConvTol = 1e-12;
% params.BarQCPConvTol = 1e-12;
params.FeasibilityTol = 1e-9;
%params.IntFeasTol = 1e-9;
% params.MarkowitzTol = 1e-4;
%params.MIPGap = 1e-12;
%params.MIPGapAbs = 1e-12;
% params.OptimalityTol = 1e-9;

% TODO: prepare inital guess for gurobi from d_lpcc and y_lpcc
if isempty(y_lpcc)
    y_lpcc = 0.0*ones(dims.n_comp,1);
end
x0 = [d_lpcc; y_lpcc];
model.start = x0;
%% solve
result = gurobi(model, params);
if isequal(result.status,'OPTIMAL')
    lpcc_solution_exists = 1;
    d_lpcc = result.x(1:dims.n_primal);
    y_lpcc = result.x(end-dims.n_comp+1:end);
    f_lin_opt = result.objval;
else
    lpcc_solution_exists = 0;
    f_lin_opt = nan;
end
lpcc_cpu_time = result.runtime;


end
