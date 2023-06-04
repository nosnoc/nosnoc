function out = extract_nonzeros_from_vector(vec)
    out = [];
    %idx = 1:length(vec);
    %out = vec(idx(logical(full(evalf(vec ~= 0)))));
    for ii=1:size(vec, 1)
        % NOTE: unfortunately this cannot be vectorized because is_zero() contains a weird bug
        x = vec(ii);
        if ~x.is_zero()
           out = vertcat(out, x);
        end
    end
end
