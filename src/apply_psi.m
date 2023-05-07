function result = apply_psi(comp_pairs, psi, sigma, scale)
    if ~exist('scale', 'var')
        scale = 1;
    end

    if isscalar(scale)
        scale = scale*ones(size(comp_pairs, 1), 1);
    elseif size(comp_pairs, 1) ~= size(scale, 1)
        error('nosnoc: scale is incorrect.')
    end
    result = [];
    for ii = 1:size(comp_pairs, 1)
        x = comp_pairs(ii, 1);
        y = comp_pairs(ii, 2);
        % NOTE: this is a bit of a hack to identify placeholder pairs
        if x.is_zero() && y.is_zero()
            result = vertcat(result, 0);
        else
            result = vertcat(result, psi(x,y,sigma*(scale(ii))));
        end
    end
end
