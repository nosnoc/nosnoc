function [ind_scalar,ind_nonscalar, ind_map] = find_nonscalar(g,w)
    % This function returns the indices of g which contain elements of the symbolic vector w, along with its
    % nonscalar indicies and the ind
    import casadi.*
    ind_g_fun = Function('ind_G', {w}, {g.jacobian(w)});
    [ind_g1,ind_g2] = find(full(sparse(DM(ind_g_fun(w) == 1) == 1)) | isnan(sparse(DM(ind_g_fun(w) == 1))));

    uniq_ind_g1 = unique(ind_g1);
    c=groupcounts(ind_g1);
    ind_scalar=uniq_ind_g1(find(c==1));
    ind_nonscalar = setdiff(1:length(g),ind_scalar)';
    ind_map = ind_g2(find(ismember(ind_scalar,ind_g1)));
end
