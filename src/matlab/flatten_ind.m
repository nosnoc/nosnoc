function flat = flatten_ind(ind)
    ind = ind.';
    flat = vertcat(ind{:});
end
