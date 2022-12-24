function inc = increment_indices(ind, len)
    inc = cellfun(@(in) in+len, ind);
end
