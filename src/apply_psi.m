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
    % TODO: for some reason casadi also wants outputs to be scalar to allow for
    % multiple argument evaluation. This is a problem for 2-sided relaxations so
    % here is a mild hack to allow it to still work.
    if (size(psi(sigma,sigma,sigma),2) == 1)
        result = psi(comp_pairs(:,1), comp_pairs(:,2), sigma*scale);
    else
        result = psi(comp_pairs(:,1)', comp_pairs(:,2)', sigma*scale);
        result = reshape(result, [length(temp),2])';
    end
    
    % Note: This is done to maintain sparsity of the outputs. For some reason
    % matrix application of functions does not preserve sparsity which
    % is either a bug or a feature I don't understand
    result = result.*sparsity;
end
