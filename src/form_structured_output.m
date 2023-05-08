function [results]= form_structured_output(problem, w_opt, name, results)
    ind = problem.(strcat('ind_', name));

    % generate structured output
    opt_s = cellfun(@(idx) w_opt(idx), ind(:,:,end), 'uni', 0);
    
    lens = cellfun(@(c) length(c), opt_s);
    if all(lens==0)
        % leave and do nothing
        return
    end

    
    len = max(lens, [], 'all');
    results.(strcat(name, '_opt_s')) = opt_s;

    % generate elementwise structured output
    i_opt = cell(len, 1);
    i_opt_flat = cell(len, 1);
    for ii = 1:len
        i_opt{ii} = cellfun(@(vec) vec(ii), opt_s(~cellfun('isempty', opt_s)));
        i_opt_flat{ii} = reshape(transpose(i_opt{ii}), prod(size(i_opt{ii})), 1);
    end

    results.(strcat(name, '_i_opt')) = i_opt;
    results.(strcat(name, '_i_opt_flat')) = i_opt_flat;

    % generate old output
    opt_ind = ind(:,:,end);
    opt_ind = opt_ind(~cellfun('isempty', opt_ind)); % hack to handle nonimpacts
    opt_ind = sort_ind_sets(opt_ind(:));
    opt_vals = cellfun(@(ind) w_opt(ind), opt_ind, 'uni', false);

    opt = horzcat(opt_vals{:});
    results.(strcat(name, '_opt')) = opt;

    opt_extended_ind = ind(:,:,:);
    opt_extended_ind = opt_extended_ind(~cellfun('isempty', opt_extended_ind)); % hack to handle differential mode
    
    opt_extended_ind = sort_ind_sets(opt_extended_ind(:));
    opt_extended_vals = cellfun(@(ind) w_opt(ind), opt_extended_ind, 'uni', false);

    opt_extended = horzcat(opt_extended_vals{:});
    results.(strcat(name, '_opt_extended')) = opt;
end
