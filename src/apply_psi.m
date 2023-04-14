function result = apply_psi(comp_pairs, psi, sigma)
    result = [];
    for ii = 1:size(comp_pairs, 1)
        x = comp_pairs(ii, 1);
        y = comp_pairs(ii, 2);
        result = vertcat(result, psi(x,y,sigma));
    end
end
