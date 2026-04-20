function [ind_scalar,ind_nonscalar, ind_map] = find_nonscalar(g,w,p)
    % Implemented by Anton Pozharskiy within nosnoc: https://github.com/nosnoc/nosnoc
    % This function returns the indices of g which contain elements of the symbolic vector w, along with its
% nonscalar indicies and the ind_map which are the indices in w which correspond to scalar elements of g.
    if ~exist('p')
         if strcmp(class(w),'casadi.SX')
            p = SX([]);
         else
            p = MX([]);
         end
    end
    import casadi.*
    % Take jacobian of g with respect to w.
    ind_g_fun = Function('ind_G', {w,p}, {g.jacobian(w)});
    if strcmp(class(w),'casadi.SX')
        p_val = SX(ones(length(p),1));
    else
        p_val = MX(ones(length(p),1));
    end
    % HERE BE DRAGONS:
    % `ind_g_fun(w) == 1` creates a symbolic array where literal `1` are ones, other constants are zeros and
    % symbolics are nans when converted to a casadi.DM.
    %
    % We want to find rows which contain only exactly one nonstructural 0, 1, or nan.
    % We get this from using find on the sparsity patern of the DM.
    if strcmp(class(w),'casadi.SX')
        sp = DM(ind_g_fun(w,p) == 1).sparsity;
    else
        % Create a function from MX expression
        % f = Function('f', {w, p}, {ind_g_fun(w, p) == 1});
        % % Evaluate with appropriate inputs (zeros work for sparsity detection)
        % w_val = DM.zeros(w.size());
        % result = f(w_val, p_val);
        % % Get sparsity
        % sp = result.sparsity();
        % Get Jacobian sparsity (this works directly with MX)
        mx_expr = ind_g_fun(w, p);
        sp = mx_expr.sparsity();
    end
    sub = sp.find;
    [ind_g1, ind_g2] = ind2sub(size(ind_g_fun(w,p)),sub);
    % transpose because groupcounts expects column vector
    [cu,ia,ic]=unique(ind_g1,"stable");
    [c,cg]=groupcounts(ind_g1');
    all_nonscalar = cg(find(c==1));
    [ind_scalar,icu,icg]=intersect(ind_g1,all_nonscalar,"stable");
    % transpose so 0x1 vectors are correctly handled externally
    ind_map = ind_g2(icu);

    % HERE BE MORE DRAGONS:
    % We also need to check for constant offsets. In future we may want to capture the offset in lbw but for now just assume nonscalar
    check = [];
    if strcmp(class(w),'casadi.SX')
        if length(ind_scalar) > 0
            check = full(DM(g(ind_scalar)-w(ind_map))) == 0;
        end
    else
        if length(ind_scalar) > 0
            diff_fun = Function('diff_check', {w}, {g(ind_scalar) - w(ind_map)});
            % You'll need numeric values here - either current values or zeros for structure check
            w_val = DM.ones(w.size());
            % Evaluate and check
            diff_result = diff_fun(w_val);
            check = full(diff_result) == 0;
        end
    end
    ind_scalar = ind_scalar(check);
    ind_map = ind_map(check);

    ind_nonscalar = setdiff(1:length(g),ind_scalar)';

    % Avoid 1x0 vs 0x1 issues
    if numel(ind_scalar) == 0
        ind_scalar = [];
    end
    if numel(ind_nonscalar) == 0
        ind_nonscalar = [];
    end
    if numel(ind_map) == 0
        ind_map = [];
    end
end
