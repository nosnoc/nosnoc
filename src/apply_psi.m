function result = apply_psi(comp_pairs, psi, sigma, scale)
    import casadi.*
    if ~exist('scale', 'var')
        scale = 1;
    end
    if isempty(comp_pairs)
        result = SX([]);
        return
    end
    temp = comp_pairs(:,1);
    sparsity = temp.sparsity();
    result = psi(comp_pairs(:,1), comp_pairs(:,2), sigma*scale);
    % Note: This is done to maintain sparsity of the outputs. For some reason
    % matrix application of functions does not preserve sparsity which
    % is either a bug or a feature I don't understand
    result = result.*sparsity;
end
