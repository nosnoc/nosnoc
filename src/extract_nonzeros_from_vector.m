function out = extract_nonzeros_from_vector(vec)
    out = [];
    for ii=1:size(vec, 1)
        % NOTE: unfortunately this cannot be vectorized because is_zero() contains a weird bug
        x = vec(ii);
        if ~x.is_zero()
           out = vertcat(out, x);
        end
    end
end
