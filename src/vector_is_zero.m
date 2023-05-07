function out = vector_is_zero(vec)
    out = [];
    for ii=1:size(vec, 1)
        % NOTE: unfortunately this cannot be vectorized because is_zero() contains a weird bug
        x = vec(ii);
        out = vertcat(out, ~(x.is_zero()));
    end
end
