function nonzero = count_nonzero(vec)
    nonzero = 0;
    for ii=1:size(vec, 1)
        % NOTE: unfortunately this cannot be vectorized because is_zero() contains a weird bug
        x = vec(ii);
        if ~x.is_zero()
           nonzero = nonzero + 1;
        end
    end
end
