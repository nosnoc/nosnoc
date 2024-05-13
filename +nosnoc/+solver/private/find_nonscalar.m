function [ind_scalar,ind_nonscalar, ind_map] = find_nonscalar(g,w)
    % This function returns the indices of g which contain elements of the symbolic vector w, along with its
    % nonscalar indicies and the ind_map which are the indices in w which correspond to scalar elements of g.
    import casadi.*
    % Take jacobian of g with respect to w.
    ind_g_fun = Function('ind_G', {w}, {g.jacobian(w)});
    % HERE BE DRAGONS:
    % `ind_g_fun(w) == 1` creates a symbolic array where literal `1` are ones, other constants are zeros and
    % symbolics are nans when converted to a casadi.DM.
    %
    % We want to find rows which contain only exactly one nonstructural 0, 1, or nan.
    % We get this from using find on the sparsity patern of the DM.
    sp = DM(ind_g_fun(w) == 1).sparsity;
    sub = sp.find;
    [ind_g1, ind_g2] = ind2sub(size(ind_g_fun(w)),sub);
    % transpose because groupcounts expects column vector
    c=groupcounts(ind_g1');
    ind_scalar=ind_g1(find(c==1));
    % transpose so 0x1 vectors are correctly handled externally
    ind_map = ind_g2(find(ismember(ind_scalar,ind_g1)));

    % HERE BE MORE DRAGONS:
    % We also need to check for constant offsets. In future we may want to capture the offset in lbw but for now just assume nonscalar
    check = full(DM(g(ind_scalar)-w(ind_map))) == 0;
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
