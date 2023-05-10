function sorted = sort_ind_sets(ind_set)
    initial_shape = size(ind_set);
    first_ind = cellfun(@(iset) iset(1), ind_set);
    
    [~,sorted_idx] = sort(first_ind);
    sorted = ind_set(sorted_idx);
end
