function [results]= form_structured_output(problem, w_opt, name, results)
    ind = problem.(strcat('ind_', name));

    % generate structured output
    opt_s = cellfun(@(idx) w_opt(idx), ind(:,:,end), 'uni', 0);
    
    lens = cellfun(@(c) length(c), opt_s);
    if all(lens==0)
        % leave and do nothing
        return
    end

    results.structured.(name) = opt_s;

    temp = opt_s';
    flat = horzcat(temp{:});
    results.(name) = flat;

    % generate extended
    opt_extended = cellfun(@(idx) w_opt(idx), sort_ind_sets(ind(:)), 'uni', 0);

    temp = opt_extended';
    results.extended.(name) = horzcat(temp{:});
end
